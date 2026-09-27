#!/usr/bin/env python3
"""
Modal analysis of the polychora as thin-walled HOLLOW 4D bodies.

Usage:
  python modal_hollow.py                    # every polychoron in topology_output/
  python modal_hollow.py pen tes hi         # selected polychora

Output: modal_hollow_output/<name>.json, same format as modal_output (see README.md).

Physics (thin-wall limit)
-------------------------
The boundary of the polychoron is a shell of wall thickness h made of its 3D cells. Each
cell is a flat 3D plate ("hyperplate") that bends in the 4th direction, along its normal,
with Kirchhoff energy (plane stress, sigma_nn = 0)

    U = (h^3/24) int_cell [lam* (Lap w)^2 + 2 mu |grad grad w|^2],  lam* = 2 lam mu / (lam + 2 mu).

A ridge (polygon shared by two cells) cannot move: any motion in its 2D normal plane would
stretch at least one of the two walls in-plane, which is infinitely stiffer than bending as
h -> 0. So every wall's deflection vanishes on all its faces. Neighbouring walls are welded:
the hinge rotation across a ridge (the slope, taken in the ridge's normal plane with the
orientation sign of n_A^m_A versus n_B^m_B) is continuous, imposed by a stiff penalty.
The bending spectrum is then exactly linear in h: frequency(h) = frequency * h / h_ref.

Discretisation (run_dist)
-------------------------
Per cell, polynomials (Legendre products, total degree <= p, local 3D coordinates) times
bubbles that vanish on every face (CellShape): phi^2 * poly with phi = (sum_f d_f^-2)^(-1/2)
(zero value and slope: the interior) plus phi * poly and b * poly, b = prod_f tanh(d_f/delta),
for the hinge rotations at the welds. Cell matrices once per shape (symmetric Gauss rule on
the flag tetrahedra, basis orthonormalised by QR), carried to every congruent cell by an
orthogonal map found from the vertex sets. Large (soft) walls: slopes 6, interior 8; small
walls: slopes 6, interior 4 if they resonate in the computed range, else phi slopes only and
interior 2. The problem is block-diagonalised by the commuting mirror symmetries (Z2)^k
(analyse_hollow_sym); blocks related by a mirror permutation are computed once. Modes up to
5 x f1 (plus the band reaching over it). Error indicator: change of two symmetry blocks with
every degree lowered by 2 (pessimistic).

Also: --isolated (one large wall + its neighbours, band means only, for prahi, prix,
gidpixhi; run_isolated) and --bubble tanh (the first version: flat-top bubble only, no
symmetry reduction, too high by up to 33 % for walls with many faces; run_tanh).

Units: circumradius R = 1, rho = 1, mu = 1, lambda = mu (as modal_output), reference wall
thickness h_ref = 0.06 R. f_Hz = frequency * (h / h_ref) * c_s / R.
"""
import argparse
import glob
import itertools
import json
import math
import os
import sys
import time

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.spatial import cKDTree

from modal_analysis import Monomials, Polytope, _unit, group_multiplets, sig

FORMAT = "polychoron-modal/1"
H_REF = 0.06


def gauss01(q):
    t, w = np.polynomial.legendre.leggauss(q)
    return (t + 1) / 2, w / 2


def tet_rule(q):
    """Collapsed Gauss rule on the ordered simplex {1 >= z1 >= z2 >= z3 >= 0} (volume 1/6)."""
    t, w = gauss01(q)
    T = np.stack(np.meshgrid(t, t, t, indexing="ij"), -1).reshape(-1, 3)
    W = np.prod(np.stack(np.meshgrid(w, w, w, indexing="ij"), -1).reshape(-1, 3), 1)
    return np.cumprod(T, axis=1), W * T[:, 0] ** 2 * T[:, 1]


def tri_rule(q):
    t, w = gauss01(q)
    U, V = np.meshgrid(t, t, indexing="ij")
    return np.stack([U * (1 - V), U * V], -1).reshape(-1, 2), (np.outer(w, w) * U).ravel()


def mirror_group(poly, tol=1e-6):
    """Commuting symmetries (Z2)^k of the polytope, as (generators, elements): a largest set
    of mutually orthogonal mirror hyperplanes (candidates: the edge directions and the
    coordinate axes), plus the central inversion if it is a symmetry and not generated.
    Element m (bit mask over the generators) is the product of the generators in m."""
    X = poly.xi
    tree = cKDTree(X)
    cellset = {frozenset(c.tolist()) for c in poly.cells}

    def invariant(T):
        dist, perm = tree.query(X @ T.T)
        return dist.max() < tol and all(frozenset(perm[c].tolist()) in cellset
                                        for c in poly.cells)

    D = X[poly.edges[:, 1]] - X[poly.edges[:, 0]]
    D = np.concatenate([D / np.linalg.norm(D, axis=1, keepdims=True), np.eye(4)])
    D *= np.where(D[np.arange(len(D)), np.abs(D).argmax(axis=1)] < 0, -1.0, 1.0)[:, None]
    cand = np.unique(np.round(D, 6), axis=0)
    normals = [n / np.linalg.norm(n) for n in cand
               if invariant(np.eye(4) - 2 * np.outer(n, n) / (n @ n))]
    best = []

    def extend(chosen, start):
        nonlocal best
        if len(chosen) > len(best):
            best = list(chosen)
        for j in range(start, len(normals)):
            if len(best) == 4:
                return
            if all(abs(normals[j] @ normals[i]) < 1e-6 for i in chosen):
                extend(chosen + [j], j + 1)

    extend([], 0)
    gens = []
    for j in best:                                      # exact to rounding: the orthogonal map
        R = np.eye(4) - 2 * np.outer(normals[j], normals[j])   # that best maps the vertices
        _, perm = tree.query(X @ R.T)                   # onto their images (Procrustes)
        u, _, vt = np.linalg.svd(X.T @ X[perm])
        gens.append((u @ vt).T)
    if len(gens) < 4 and invariant(-np.eye(4)):
        gens.append(-np.eye(4))
    elements = []
    for m in range(2 ** len(gens)):
        g = np.eye(4)
        for j in range(len(gens)):
            if m >> j & 1:
                g = g @ gens[j]
        elements.append(g)
    return gens, elements


def mirror_permutations(poly, gens, tol=1e-6):
    """Symmetries g that permute the mirror generators by conjugation, g R_i g^-1 = R_pi(i)
    (the central inversion, if a generator, stays): list of (pi, g). They map the block of
    character chi onto the block of chi o pi^-1, which then has the same spectrum."""
    X = poly.xi
    tree = cKDTree(X)
    cellset = {frozenset(c.tolist()) for c in poly.cells}
    mir = [j for j, g in enumerate(gens) if abs(np.trace(g) - 2) < 1e-9]
    normals = {j: np.linalg.eigh(gens[j])[1][:, 0] for j in mir}   # eigenvalue -1 first
    out = []
    for perm in itertools.permutations(mir):
        pi = {j: j for j in range(len(gens))}
        pi.update(dict(zip(mir, perm)))
        for signs in itertools.product([1.0, -1.0], repeat=len(mir)):
            g = sum(s * np.outer(normals[pi[j]], normals[j]) for j, s in zip(mir, signs))
            if len(mir) < 4:                            # complete the map on the rest
                Nm = np.array([normals[j] for j in mir])
                rest = sla.null_space(Nm)
                g = g + rest @ rest.T
            dist, pm = tree.query(X @ g.T)
            if dist.max() < tol and all(frozenset(pm[c].tolist()) in cellset for c in poly.cells):
                u, _, vt = np.linalg.svd(X.T @ X[pm])
                out.append((pi, (u @ vt).T))
                break
    return out


def face_slopes(poly, f, cells, shapes, cell_labels, frames, Qs, U, W):
    """Quadrature weights on ridge f and the hinge slopes of the given cells (one or both of
    the ridge's cells) at its points; with both, the second carries the orientation sign, so
    that the weld energy is kr * int (g_1 a_1 + g_2 a_2)^2."""
    verts, xf = poly.xi[poly.faces[f]], poly.face_ctr[f]
    _, _, vt = np.linalg.svd(verts - xf)
    ang = np.arctan2((verts - xf) @ vt[1], (verts - xf) @ vt[0])
    cyc = verts[np.argsort(ang)]
    X, Wq = [], []
    for k in range(len(cyc)):
        e1, e2 = cyc[k] - xf, cyc[(k + 1) % len(cyc)] - xf
        area2 = math.sqrt(max(e1 @ e1 * (e2 @ e2) - (e1 @ e2) ** 2, 0))
        X.append(xf + U[:, :1] * e1 + U[:, 1:] * e2)
        Wq.append(W * area2)
    X, Wq = np.concatenate(X), np.concatenate(Wq)
    gs, nm = [], []
    for Xc in cells:
        sh = shapes[cell_labels[Xc]]
        o, T, r = frames[Xc]
        s = ((X - o) @ T.T / r) @ Qs[Xc]               # reference-cell coordinates
        # the same face in the reference cell: the reference face whose plane contains s
        mref = s @ sh.m.T - sh.d[None, :]
        kf = int(np.argmin(np.abs(mref).max(axis=0)))
        gs.append(sh.face_slope(kf, s) / r)
        m_in = -(Qs[Xc] @ sh.m[kf])                      # inward, in the cell's own local frame
        nm.append((poly.cell_normal[Xc], T.T @ m_in))
    if len(cells) == 2:
        (nA, mA), (nB, mB) = nm
        gs[1] = -np.sign(nA @ nB * (mA @ mB) - nA @ mB * (mA @ nB)) * gs[1]
    return Wq, gs


def cell_frame(poly, c):
    """Centroid o, orthonormal tangent basis T (3 x 4, orthogonal to the normal), circumradius."""
    P = poly.xi[poly.cells[c]]
    o = P.mean(axis=0)
    n = poly.cell_normal[c]
    _, _, vt = np.linalg.svd(P - o)
    T = vt[:3] - np.outer(vt[:3] @ n, n)
    T = np.linalg.qr(T.T)[0].T
    return o, T, float(np.linalg.norm(P - o, axis=1).max())


class CellShape:
    """Reference cell of one shape: face planes, bubble basis, bending/mass matrices.
    Local coordinates s = T (x - o) / r (r = cell circumradius); matrices are per unit
    thickness factors: K = kb * h^3 / 12 / r, M = G * h * r^3 with the returned kb, G.

    Two bubbles: the flat-top b = prod_f tanh(d_f / delta) (smooth, vanishes linearly on each
    face, its slope dies out towards the face's edges), and phi = (sum_f d_f^-2)^(-1/2), a
    smooth distance to the face planes (phi ~ d_f at face f, no layer width).
    bubble "dist" (default): phi^2 * poly(degree p), zero value and slope on every face (the
    interior; a clamped wall converges to < 1 % at degree 6), plus, unless clamped,
    phi * poly(degree p_slope), which carries the hinge rotations at the faces. Near a cell
    edge phi * poly behaves like r, the exact solution like r^(pi/alpha): best for the obtuse
    edges of Archimedean cells (simply supported tetrahedron: +0.14 % at p_slope = 6).
    bubble "both": the rotations carried by b * poly and phi * poly (slightly better per
    degree, more unknowns). bubble "tanh": b * poly(degree p) only (first version; a clamped
    or stiffly welded wall needs a poly vanishing on all faces, which converges slowly).
    The basis is orthonormalised (G = I) by a QR of the weighted values."""

    def __init__(self, poly, c, p, lam_s, mu, q=None, delta_frac=None, bubble="dist",
                 clamped=False, p_slope=None):
        self.o, self.T, self.r = cell_frame(poly, c)
        self.bubble = bubble
        ps = p if p_slope is None else p_slope
        if bubble == "tanh":
            self.groups = [("tanh", 1, p)]                 # (bubble, power, degree)
        elif clamped:
            self.groups = [("dist", 2, p)]
        elif bubble == "both":
            self.groups = [("tanh", 1, ps), ("dist", 1, ps), ("dist", 2, p)]
        else:
            self.groups = [("dist", 1, ps), ("dist", 2, p)]
        self.p = max(pg for _, _, pg in self.groups)
        # polynomials: Legendre products P_a1(s1) P_a2(s2) P_a3(s3), total degree <= p (the same
        # space as the monomials, better conditioned on the cell, |s| <= 1)
        self.E = Monomials(3, self.p).upto(self.p)
        deg = self.E.sum(axis=1)
        self.cols = [np.flatnonzero(deg <= pg) for _, _, pg in self.groups]
        loc = self.local
        self.V = loc(poly.xi[poly.cells[c]])
        self.faces = []                                # (outward unit normal, distance) local
        for f in poly.cell_faces[c]:
            fv = poly.xi[poly.faces[f]]
            fb = np.linalg.svd(fv - fv.mean(0))[2][:2]
            v = (poly.face_ctr[f] - self.o)
            v = v - fb.T @ (fb @ v)
            self.faces.append((self.T @ _unit(v), np.linalg.norm(v) / self.r))
        self.m = np.array([m for m, _ in self.faces])
        self.d = np.array([d for _, d in self.faces])
        if delta_frac is None:                         # auto: nearly polynomial bubble for few
            F = len(self.d)                            # faces (fastest convergence), flat-top
            delta_frac = 3.0 if F <= 8 else 1.0 if F <= 12 else 0.5   # for many
        self.delta = self.d.min() * delta_frac
        self.norm = {"tanh": np.prod(np.tanh(self.d / self.delta)), "dist": 1.0}
        self.norm["dist"] = float(self._bubble(np.zeros((1, 3)), "dist")[0][0])
        # quadrature on the flag tetrahedra (centroid, face centre, edge endpoints); the dist
        # basis uses both orders of the edge endpoints, so that the rule has the cell's full
        # symmetry (the collapsed Gauss rule is not symmetric in them)
        qq = q or self.p + 6
        Z, W = tet_rule(qq)
        S, Wt = [], []
        for f in poly.cell_faces[c]:
            for e in poly.face_edges[f]:
                a, bv = poly.xi[poly.edges[e]]
                for ends in ((a, bv), (bv, a)) if bubble != "tanh" else ((a, bv),):
                    Vt = loc(np.stack([self.o, poly.face_ctr[f], *ends]))
                    D = Vt[1:] - Vt[:-1]
                    S.append(Vt[0] + Z @ D)
                    Wt.append(W * abs(np.linalg.det(D)) / (2 if bubble != "tanh" else 1))
        S, Wt = np.concatenate(S), np.concatenate(Wt)
        nb = sum(len(cl) for cl in self.cols)
        chunk = max(500, 2_000_000 // nb)
        # orthonormal basis without forming the (ill-conditioned) Gram matrix: R from a QR of
        # the weighted values, accumulated over chunks (TSQR); basis psi = raw @ R^-1. The
        # error then grows with sqrt(cond G) only, so the symmetries stay exact.
        R = np.zeros((0, nb))
        for s0 in range(0, len(S), chunk):
            phi = self.eval(S[s0:s0 + chunk], values_only=True)
            R = np.linalg.qr(np.vstack([R, np.sqrt(Wt[s0:s0 + chunk])[:, None] * phi]), mode="r")
        R *= np.sign(np.diag(R))[:, None]
        self.R = R
        kb = np.zeros((nb, nb))
        for s0 in range(0, len(S), chunk):
            Ws = Wt[s0:s0 + chunk]
            _, _, hess = self.eval(S[s0:s0 + chunk])
            ho = {}
            for i in range(3):
                for j in range(i, 3):
                    ho[i, j] = sla.solve_triangular(R, hess[i, j].T, trans="T").T  # hess @ R^-1
            lap = ho[0, 0] + ho[1, 1] + ho[2, 2]
            kb += lam_s * lap.T @ (Ws[:, None] * lap)
            for (i, j), hij in ho.items():
                kb += (2 if i == j else 4) * mu * hij.T @ (Ws[:, None] * hij)
        self.volume = float(Wt.sum())
        self.kb, self.G = (kb + kb.T) / 2, np.eye(nb)
        self.size = nb
        self.samples = S[np.random.default_rng(0).choice(len(S), min(len(S), 4 * nb), replace=False)]

    def orth(self, A):
        """Raw-basis values (rows) -> values of the orthonormal basis: A @ R^-1."""
        return sla.solve_triangular(self.R, A.T, trans="T").T

    def local(self, X):
        return (X - self.o) @ self.T.T / self.r

    def rep_matrix(self, B):
        """Orthogonal matrix D acting on this shape's (orthonormal) coefficients as the map
        w -> w o B^-1, for an orthogonal B that maps the reference cell onto itself (the
        bubble is invariant, so the space is: fitted exactly by least squares at points of
        the cell)."""
        S = self.samples
        A0 = self.orth(self.eval(S, values_only=True))
        A1 = self.orth(self.eval(S @ B, values_only=True))      # at B^-1 s
        return np.linalg.lstsq(A0, A1, rcond=None)[0]

    def _poly(self, s, derivs=False):
        """Legendre products at points s: values (m, n); with derivs also the gradient (list
        of 3) and the hessian (dict (i, j), i <= j)."""
        n = self.p + 1
        P = np.zeros((3, len(s), n))
        dP = np.zeros_like(P)
        d2P = np.zeros_like(P)
        x = s.T
        P[:, :, 0] = 1.0
        if n > 1:
            P[:, :, 1] = x
            dP[:, :, 1] = 1.0
        for k in range(1, n - 1):
            P[:, :, k + 1] = ((2 * k + 1) * x * P[:, :, k] - k * P[:, :, k - 1]) / (k + 1)
            dP[:, :, k + 1] = dP[:, :, k - 1] + (2 * k + 1) * P[:, :, k]
            d2P[:, :, k + 1] = d2P[:, :, k - 1] + (2 * k + 1) * dP[:, :, k]
        a = self.E
        F = [P[i][:, a[:, i]] for i in range(3)]
        val = F[0] * F[1] * F[2]
        if not derivs:
            return val, None, None
        D1 = [dP[i][:, a[:, i]] for i in range(3)]
        D2 = [d2P[i][:, a[:, i]] for i in range(3)]
        grad = [D1[0] * F[1] * F[2], F[0] * D1[1] * F[2], F[0] * F[1] * D1[2]]
        hess = {(0, 0): D2[0] * F[1] * F[2], (1, 1): F[0] * D2[1] * F[2],
                (2, 2): F[0] * F[1] * D2[2], (0, 1): D1[0] * D1[1] * F[2],
                (0, 2): D1[0] * F[1] * D1[2], (1, 2): F[0] * D1[1] * D1[2]}
        return val, grad, hess

    def _bubble(self, s, kind):
        """Bubble (power 1) at points s, with grad b = -b gv and the hessian of log b."""
        if kind == "tanh":
            L = self.d[None, :] - s @ self.m.T
            t = np.tanh(L / self.delta)
            sech2 = 1 - t * t
            g = sech2 / (self.delta * t)                # t'/t
            c2 = -2 * sech2 / self.delta ** 2 - g * g   # t''/t - (t'/t)^2
            return (np.prod(t, axis=1) / self.norm["tanh"], g @ self.m,
                    np.einsum("mf,fi,fj->ijm", c2, self.m, self.m))
        u = self.d[None, :] - s @ self.m.T              # distances to the face planes
        u2 = u ** -2.0
        S = u2.sum(axis=1)
        a = u2 / S[:, None]
        v = (a / u) @ self.m
        H = (-3 * np.einsum("mf,fi,fj->ijm", a / u ** 2, self.m, self.m)
             + 2 * np.einsum("mi,mj->ijm", v, v))
        return S ** -0.5 / self.norm["dist"], v, H

    def eval(self, s, values_only=False):
        """phi (m, nb), grad (3, m, nb), hessian (3, 3, m, nb) of the raw (not orthonormalised)
        basis at interior points s (local); values_only: phi only."""
        bub = {k: self._bubble(s, k) for k in {kind for kind, _, _ in self.groups}}
        if values_only:
            mono = self._poly(s)[0]
            return np.concatenate([(bub[k][0] ** e)[:, None] * mono[:, cols]
                                   for (k, e, _), cols in zip(self.groups, self.cols)], axis=1)
        mono, dmono, d2mono = self._poly(s, derivs=True)
        phis, grads, hesss = [], [], []
        for (kind, e, _), cols in zip(self.groups, self.cols):
            b1, gv1, H1 = bub[kind]
            b = b1 ** e
            gv = e * gv1                                 # grad b^e = -b^e gv
            gb = -b[:, None] * gv
            Hb = b[None, None, :] * (np.einsum("mi,mj->ijm", gv, gv) + e * H1)
            mo, dm = mono[:, cols], [d[:, cols] for d in dmono]
            phis.append(b[:, None] * mo)
            grads.append(np.stack([gb[:, i:i + 1] * mo + b[:, None] * dm[i] for i in range(3)]))
            hess = np.empty((3, 3, len(s), len(cols)))
            for i in range(3):
                for j in range(i, 3):
                    hess[i, j] = (Hb[i, j][:, None] * mo + gb[:, i:i + 1] * dm[j]
                                  + gb[:, j:j + 1] * dm[i] + b[:, None] * d2mono[i, j][:, cols])
                    hess[j, i] = hess[i, j]
            hesss.append(hess)
        return (np.concatenate(phis, axis=1), np.concatenate(grads, axis=2),
                np.concatenate(hesss, axis=3))

    def clamped_omegas(self, poly, c, beta, h, mu=1.0):
        """Frequencies of this wall alone with all its faces clamped by the same slope penalty
        as the welds (the stiff limit of the welded joint), for the convergence indicator."""
        U, W = tri_rule(self.p + 5)
        P = np.zeros_like(self.kb)
        for k, f in enumerate(poly.cell_faces[c]):
            verts, xf = poly.xi[poly.faces[f]], poly.face_ctr[f]
            _, _, vt = np.linalg.svd(verts - xf)
            cyc = verts[np.argsort(np.arctan2((verts - xf) @ vt[1], (verts - xf) @ vt[0]))]
            for j in range(len(cyc)):
                e1, e2 = cyc[j] - xf, cyc[(j + 1) % len(cyc)] - xf
                area2 = math.sqrt(max(e1 @ e1 * (e2 @ e2) - (e1 @ e2) ** 2, 0)) / self.r ** 2
                g = self.face_slope(k, self.local(xf + U[:, :1] * e1 + U[:, 1:] * e2))
                P += g.T @ ((W * area2)[:, None] * g)
        K = h ** 3 / 12 / self.r * self.kb + beta * mu * h ** 3 / self.r * P
        return np.sqrt(sla.eigh(K, h * self.r ** 3 * self.G, eigvals_only=True))

    def face_slope(self, k, s):
        """Inward normal derivative (local units) of all basis functions at points s on face k."""
        mono = self._poly(s)[0]
        g = np.zeros((len(s), self.size))                # squared bubbles: zero slope
        off = 0
        for (kind, e, _), cols in zip(self.groups, self.cols):
            if e == 1 and kind == "tanh":
                L = self.d[None, :] - s @ self.m.T
                L[:, k] = 0.0
                t = np.tanh(L / self.delta)
                t[:, k] = 1.0
                db = np.prod(t, axis=1) / (self.delta * self.norm["tanh"])
                g[:, off:off + len(cols)] = db[:, None] * mono[:, cols]
            elif e == 1:
                # phi = d_k (1 + sum_f (d_k/d_f)^2)^(-1/2): slope exactly 1 on the open face k
                g[:, off:off + len(cols)] = mono[:, cols] / self.norm["dist"]
            off += len(cols)
        return self.orth(g)


def congruence(ref, V, adjacency_ref, adjacency, tol=1e-6):
    """Orthogonal Q (3x3) with Q @ ref.V[i] = V[perm[i]] for all vertices (congruent cells)."""
    tree = cKDTree(V)
    a = 0
    nb = sorted(adjacency_ref[a])
    for b_ in nb:
        for c_ in nb:
            if c_ == b_:
                continue
            A = np.stack([ref.V[a], ref.V[b_], ref.V[c_]]).T
            if abs(np.linalg.det(A)) < 1e-6:
                continue
            Ainv = np.linalg.inv(A)
            dbc = np.linalg.norm(ref.V[b_] - ref.V[c_])
            for a2 in range(len(V)):
                for b2 in adjacency[a2]:
                    for c2 in adjacency[a2]:
                        if c2 == b2 or abs(np.linalg.norm(V[b2] - V[c2]) - dbc) > tol:
                            continue
                        Q = np.stack([V[a2], V[b2], V[c2]]).T @ Ainv
                        if np.abs(Q.T @ Q - np.eye(3)).max() > 1e-6:
                            continue
                        dist, _ = tree.query(ref.V @ Q.T)
                        if dist.max() < tol:
                            return Q
            return None
    return None


def analyse_hollow(poly, cell_labels, degree, beta, n_modes, lam=1.0, mu=1.0, bubble="dist",
                   delta_frac=None, cells=None, clamp_outside=None, shapes=None, verbose=True):
    """Welded thin-walled hollow polychoron. Returns frequencies (h = H_REF), multiplets,
    per-cell ∫ w^2 contributions for gains.
    degree: per cell label, p or (p_slope, p) (see CellShape); shapes: prebuilt CellShapes.
    cells: compute only these cells (default all); a ridge to a cell outside is clamped if
    clamp_outside(cell) is true, else hinged (free rotation)."""
    t0 = time.time()
    lam_s = 2 * lam * mu / (lam + 2 * mu)
    h = H_REF
    cells = list(range(len(poly.cells))) if cells is None else list(cells)
    where = {c: i for i, c in enumerate(cells)}
    labels = sorted({cell_labels[c] for c in cells})
    shapes = dict(shapes or {})
    for lab in labels:
        if lab in shapes:
            continue
        rep = cell_labels.index(lab)
        dg = degree[lab]
        p, ps = (dg[1], dg[0]) if isinstance(dg, tuple) else (dg, None)
        shapes[lab] = CellShape(poly, rep, p, lam_s, mu, delta_frac=delta_frac, bubble=bubble,
                                p_slope=ps)
    # cell adjacency (vertex graph inside each cell) for congruence
    def cell_adj(c):
        verts = list(poly.cells[c])
        idx = {v: i for i, v in enumerate(verts)}
        adj = [set() for _ in verts]
        for f in poly.cell_faces[c]:
            for e in poly.face_edges[f]:
                u, v = poly.edges[e]
                adj[idx[u]].add(idx[v]); adj[idx[v]].add(idx[u])
        return adj
    Qs, frames_ = {}, {}
    ref_adj = {lab: cell_adj(cell_labels.index(lab)) for lab in labels}
    for c in cells:
        sh = shapes[cell_labels[c]]
        o, T, r = cell_frame(poly, c)
        V = (poly.xi[poly.cells[c]] - o) @ T.T / r
        Q = congruence(sh, V, ref_adj[cell_labels[c]], cell_adj(c))
        if Q is None:
            raise RuntimeError(f"{poly.name}: cell {c} not congruent to its reference")
        Qs[c] = Q
        frames_[c] = (o, T, r)
    t1 = time.time()

    sizes = [shapes[cell_labels[c]].size for c in cells]
    offs = np.concatenate([[0], np.cumsum(sizes)])
    N = int(offs[-1])
    Kb, Mb = {}, {}
    for i, c in enumerate(cells):
        sh = shapes[cell_labels[c]]
        Kb[i, i] = h ** 3 / 12 / sh.r * sh.kb
        Mb[i, i] = h * sh.r ** 3 * sh.G
    ell = float(np.mean([shapes[l].r for l in labels]))
    kr = beta * mu * h ** 3 / ell

    face_cells = [[] for _ in poly.faces]
    for c, fs in enumerate(poly.cell_faces):
        for f in fs:
            face_cells[f].append(c)
    U, W = tri_rule(max(sh.p for sh in shapes.values()) + 5)
    for f, (A, B) in enumerate(face_cells):
        inside = [X_ for X_ in (A, B) if X_ in where]
        if not inside:
            continue
        if len(inside) == 1 and not (clamp_outside is None or clamp_outside(B if A in where else A)):
            continue                                   # hinged to a cell outside
        Wq, gs = face_slopes(poly, f, inside, shapes, cell_labels, frames_, Qs, U, W)
        pairs = [(where[X_], g_) for X_, g_ in zip(inside, gs)]   # one cell: clamped
        for I, gI in pairs:
            for J, gJ in pairs:
                blk = kr * gI.T @ (Wq[:, None] * gJ)
                Kb[I, J] = Kb[I, J] + blk if (I, J) in Kb else blk
    t2 = time.time()

    def assemble(blocks):
        rows, cols, vals = [], [], []
        for (I, J), B in blocks.items():
            ii, jj = np.nonzero(B)
            rows.append((offs[I] + ii).astype(np.int32)); cols.append((offs[J] + jj).astype(np.int32))
            vals.append(B[ii, jj])
        return sp.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                             shape=(N, N))
    K, M = assemble(Kb), assemble(Mb)
    del Kb
    K = (K + K.T) / 2
    k = min(n_modes, N - 2)
    if N <= 5000:
        w, X = sla.eigh(K.toarray(), M.toarray(), subset_by_index=[0, k - 1])
    else:
        # shift-invert about 0 (K is positive definite) with a symmetric fill-reducing
        # ordering: the default COLAMD ordering of eigsh fills in badly here (GBs)
        lu = spla.splu(K.tocsc(), permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.0,
                       options=dict(SymmetricMode=True))
        OPinv = spla.LinearOperator((N, N), matvec=lu.solve, dtype=float)
        w, X = spla.eigsh(K, k=k, M=M, sigma=0.0, which="LM", OPinv=OPinv)
        o_ = np.argsort(w); w, X = w[o_], X[:, o_]
    t3 = time.time()
    omega = np.sqrt(np.maximum(w, 0))
    # int_cell w^2 per cell and mode (mass matrix blocks), for the gains
    wsq = np.zeros((len(cells), len(omega)))
    for i in range(len(cells)):
        Xc = X[offs[i]:offs[i + 1]]
        wsq[i] = np.einsum("ak,ak->k", Xc, Mb[i, i] @ Xc) / h
    if verbose:
        print(f"  {poly.name:9s} cells={len(cells):5d} DOFs={N:6d} degrees="
              + ",".join(f"{lab.split()[-1][:4]}:{degree[lab]}" for lab in labels)
              + f"  f1={omega[0] / (2 * np.pi):.5f}  "
              f"[{t1 - t0:.0f}s shapes, {t2 - t1:.0f}s coupling, {t3 - t2:.0f}s eigen]", flush=True)
    return omega, wsq, degree, shapes


def analyse_hollow_sym(poly, cell_labels, degree, beta, omega_max, max_modes, lam=1.0, mu=1.0,
                       characters=None, shapes=None, density=0.0, verbose=True):
    """Welded hollow polychoron, block-diagonalised by the commuting symmetries (Z2)^k of
    mirror_group: one real block per character chi, spanned by chi-symmetric combinations
    over each cell orbit. A block is assembled from the rows of one representative cell per
    orbit only (K is invariant), so the full matrix is never formed.
    Computes all modes with omega <= omega_max, at most ~max_modes in total (then fewer:
    'complete' is the highest omega below which every block is complete).
    Returns omega (sorted), wsq (cells x modes, int_cell w^2 per cell), complete, shapes,
    chi index per mode, group order."""
    t0 = time.time()
    lam_s = 2 * lam * mu / (lam + 2 * mu)
    h = H_REF
    C = len(poly.cells)
    gens, G = mirror_group(poly)
    nG = len(G)
    ctr = np.array([poly.xi[c].mean(axis=0) for c in poly.cells])
    tree = cKDTree(ctr)
    perm = np.empty((nG, C), dtype=int)
    for i, g in enumerate(G):
        dist, perm[i] = tree.query(ctr @ g.T)
        if dist.max() > 1e-8:
            raise RuntimeError(f"{poly.name}: symmetry {i} does not map cells onto cells")
    rep_of, via, reps = np.full(C, -1), np.zeros(C, dtype=int), []
    for c in range(C):
        if rep_of[c] < 0:
            reps.append(c)
            for i in range(nG):
                if rep_of[perm[i, c]] < 0:
                    rep_of[perm[i, c]], via[perm[i, c]] = c, i
    size = {r: int((rep_of == r).sum()) for r in reps}
    stab = {r: [i for i in range(nG) if perm[i, r] == r] for r in reps}

    labels = sorted(set(cell_labels))
    shapes = dict(shapes or {})
    for lab in labels:
        if lab not in shapes:
            dg = degree[lab]                            # p, (p_slope, p) or (bubble, p_slope, p)
            dg = dg if isinstance(dg, tuple) else (dg, dg)
            dg = dg if len(dg) == 3 else ("dist",) + dg
            shapes[lab] = CellShape(poly, cell_labels.index(lab), dg[2], lam_s, mu,
                                    p_slope=dg[1], bubble=dg[0])

    def cell_adj(c):
        verts = list(poly.cells[c])
        idx = {v: i for i, v in enumerate(verts)}
        adj = [set() for _ in verts]
        for f in poly.cell_faces[c]:
            for e in poly.face_edges[f]:
                u, v = poly.edges[e]
                adj[idx[u]].add(idx[v]); adj[idx[v]].add(idx[u])
        return adj
    ref_adj = {lab: cell_adj(cell_labels.index(lab)) for lab in labels}
    face_cells = [[] for _ in poly.faces]
    for c, fs in enumerate(poly.cell_faces):
        for f in fs:
            face_cells[f].append(c)
    need = set(reps)
    for r in reps:
        need.update(x for f in poly.cell_faces[r] for x in face_cells[f])
    frames_, Qs = {}, {}
    for c in need:
        sh = shapes[cell_labels[c]]
        o, T, r_ = cell_frame(poly, c)
        Q = congruence(sh, (poly.xi[poly.cells[c]] - o) @ T.T / r_, ref_adj[cell_labels[c]], cell_adj(c))
        if Q is None:
            raise RuntimeError(f"{poly.name}: cell {c} not congruent to its reference")
        frames_[c], Qs[c] = (o, T, r_), Q

    Dcache = {}

    def Dmat(i, c):
        """Coefficient map of symmetry i from cell c to cell perm[i, c]."""
        d = perm[i, c]
        B = Qs[d].T @ frames_[d][1] @ G[i] @ frames_[c][1].T @ Qs[c]
        key = (cell_labels[c], tuple(np.round(B, 5).ravel()))
        if key not in Dcache:
            sh = shapes[cell_labels[c]]
            dist, pv = cKDTree(sh.V).query(sh.V @ B.T)
            if dist.max() > 1e-5:
                raise RuntimeError(f"{poly.name}: cell map {i}, {c} is no symmetry of its shape")
            u, _, vt = np.linalg.svd(sh.V.T @ sh.V[pv])   # exact symmetry of the reference
            Dcache[key] = sh.rep_matrix((u @ vt).T)
        return Dcache[key]

    # rows of K for the representative cells
    ell = float(np.mean([shapes[l].r for l in labels]))
    kr = beta * mu * h ** 3 / ell
    U, W = tri_rule(max(sh.p for sh in shapes.values()) + 5)
    rows = {r: {r: h ** 3 / 12 / shapes[cell_labels[r]].r * shapes[cell_labels[r]].kb} for r in reps}
    for r in reps:
        for f in poly.cell_faces[r]:
            A, B = face_cells[f]
            Wq, (gA, gB) = face_slopes(poly, f, [A, B], shapes, cell_labels, frames_, Qs, U, W)
            gs, go, other = (gA, gB, B) if A == r else (gB, gA, A)
            rows[r][r] = rows[r][r] + kr * gs.T @ (Wq[:, None] * gs)
            rows[r][other] = rows[r].get(other, 0) + kr * gs.T @ (Wq[:, None] * go)
    t1 = time.time()

    # characters of (Z2)^k: chi_s(element m) = (-1)^popcount(s & m). A symmetry permuting
    # the mirrors maps block s onto block pi(s) with the same spectrum: one block per orbit
    copies = {}
    if characters is None:
        chis, seen = [], set()
        perms = mirror_permutations(poly, gens)
        for s in range(nG):
            if s in seen:
                continue
            chis.append(s); seen.add(s); copies[s] = []
            for pi, g in perms:
                s2 = sum(1 << pi[j] for j in range(len(gens)) if s >> j & 1)
                if s2 not in seen:
                    seen.add(s2)
                    copies[s].append((s2, tree.query(ctr @ g.T)[1]))
    else:
        chis = characters
    om_all, wsq_all, chi_of, complete = [], [], [], np.inf
    N_total = sum(size[r] * shapes[cell_labels[r]].size for r in reps)
    # modes per unknown: the caller's estimate, then the largest seen in a block
    cell_orbit = [reps.index(rep_of[c]) for c in range(C)]
    for s in chis:
        chi = np.array([(-1) ** bin(s & m).count("1") for m in range(nG)], dtype=float)
        V = {}
        for r in reps:
            Pj = sum(chi[i] * Dmat(i, r) for i in stab[r]) / len(stab[r])
            ew, ev = np.linalg.eigh((Pj + Pj.T) / 2)
            V[r] = ev[:, ew > 0.5]
        dims = [V[r].shape[1] for r in reps]
        offs = np.concatenate([[0], np.cumsum(dims)]).astype(int)
        at = {r: k for k, r in enumerate(reps)}
        N = int(offs[-1])
        if N == 0:
            continue
        rr, cc, vv = [], [], []
        for r in reps:
            if not dims[at[r]]:
                continue
            for c2, Kb in rows[r].items():
                b = rep_of[c2]
                if not dims[at[b]]:
                    continue
                blk = (math.sqrt(size[r] / size[b]) * chi[via[c2]]) * (V[r].T @ Kb @ Dmat(via[c2], b) @ V[b])
                ii, jj = np.nonzero(np.abs(blk) > 0)
                rr.append(offs[at[r]] + ii); cc.append(offs[at[b]] + jj); vv.append(blk[ii, jj])
        K = sp.csc_matrix((np.concatenate(vv), (np.concatenate(rr), np.concatenate(cc))), shape=(N, N))
        asym = abs(K - K.T).max() / abs(K).max()
        if asym > 1e-5:                                # quadrature of phi: ~1e-8
            raise RuntimeError(f"{poly.name}: block {s} not symmetric ({asym:.1e})")
        mdiag = np.concatenate([np.full(dims[at[r]], h * shapes[cell_labels[r]].r ** 3) for r in reps])
        dinv = 1 / np.sqrt(mdiag)
        Kt = sp.csc_matrix(sp.diags(dinv) @ ((K + K.T) / 2) @ sp.diags(dinv))
        cap = int(math.ceil(max_modes * N / N_total)) + 20
        if N <= 2500:
            w, X = sla.eigh(Kt.toarray(), subset_by_value=(-np.inf, omega_max ** 2))
            if len(w) > cap:
                w, X = w[:cap], X[:, :cap]
                complete = min(complete, math.sqrt(w[-1]))
        else:
            lu = spla.splu(Kt, permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.0,
                           options=dict(SymmetricMode=True))
            OPinv = spla.LinearOperator((N, N), matvec=lu.solve, dtype=float)
            k = min(max(40, int(1.2 * density * N) + 20 if density else int(0.3 * cap)), cap, N - 2)
            while True:
                w, X = spla.eigsh(Kt, k=k, sigma=0.0, which="LM", OPinv=OPinv)
                o_ = np.argsort(w); w, X = w[o_], X[:, o_]
                if w[-1] > omega_max ** 2 or k >= min(cap, N - 2):
                    break
                k = min(int(k * 1.7) + 10, cap, N - 2)
            keep = w <= omega_max ** 2
            if not keep.all():
                w, X = w[keep], X[:, keep]
            else:
                complete = min(complete, math.sqrt(w[-1]))
            del lu
        density = max(density, len(w) / N)
        Y = X * dinv[:, None]                           # M-normalised coefficients
        wo = np.empty((len(reps), len(w)))
        for k_, r in enumerate(reps):
            y = Y[offs[k_]:offs[k_ + 1]]
            wo[k_] = shapes[cell_labels[r]].r ** 3 * np.einsum("ak,ak->k", y, y) / size[r]
        om = np.sqrt(np.maximum(w, 0))
        wcell = wo[cell_orbit]                          # orbit value for every cell
        om_all.append(om); wsq_all.append(wcell); chi_of.append(np.full(len(w), s))
        for s2, pg in copies.get(s, []):                # same spectrum, cells permuted
            w2 = np.empty_like(wcell)
            w2[pg] = wcell
            om_all.append(om); wsq_all.append(w2); chi_of.append(np.full(len(w), s2))
        if verbose:
            print(f"      block {s:2d}/{nG} (x{1 + len(copies.get(s, []))}): {N:6d} unknowns, "
                  f"{len(w):5d} modes  [{time.time() - t0:.0f}s]", flush=True)
    omega = np.concatenate(om_all)
    o_ = np.argsort(omega, kind="stable")
    omega = omega[o_]
    wsq = np.concatenate(wsq_all, axis=1)[:, o_]
    t2 = time.time()
    if verbose:
        print(f"  {poly.name:9s} cells={C:5d} group (Z2)^{len(gens)}, {len(reps)} orbits, "
              f"{N_total} unknowns, degrees="
              + ",".join(f"{lab.split()[-1][:4]}:{degree[lab]}" for lab in labels)
              + f"  f1={omega[0] / (2 * np.pi):.5f}  [{t1 - t0:.0f}s setup, {t2 - t1:.0f}s blocks]",
              flush=True)
    return omega, wsq, complete, shapes, np.concatenate(chi_of)[o_], nG


def complete_bands(freq, top=0.4):
    """Number of leading multiplets that form complete bands. A hollow spectrum is banded
    (one mode per wall and wall mode), and the computed set or the frequency limit may end
    inside a band, whose modes would then be under-represented. Drop the last band: cut at
    the last relative gap in the top 40 % of the range that is at least half the widest gap
    there (and >= 3 %) -- a band boundary, not a split inside a coupled band."""
    f = np.asarray(freq)
    if len(f) < 3:
        return len(f)
    gaps = (f[1:] - f[:-1]) / f[:-1]
    ok = f[:-1] >= (1 - top) * f[-1]
    if not ok.any():
        return len(f)
    big = np.flatnonzero(ok & (gaps >= max(0.03, 0.5 * gaps[ok].max())))
    return int(big[-1]) + 1 if len(big) else len(f)


def soft_shapes(poly, cell_labels, stiff_ratio=3.0):
    """Cell shapes whose walls are soft (in the audible band): inradius within
    sqrt(stiff_ratio) of the largest (bending frequency ~ 1 / size^2)."""
    inr = {}
    for lab in set(cell_labels):
        c = cell_labels.index(lab)
        o = poly.xi[poly.cells[c]].mean(axis=0)
        inr[lab] = min(np.linalg.norm(poly.face_ctr[f] - o) for f in poly.cell_faces[c])
    big = max(inr.values())
    return {lab for lab in inr if (big / inr[lab]) ** 2 < stiff_ratio}


def choose_degrees(poly, cell_labels, max_dofs, stiff_ratio=3.0):
    """Highest degree for the soft cells (and 2 or 1 for stiff ones) within the DOF budget;
    also returns how many walls are soft."""
    counts = {lab: cell_labels.count(lab) for lab in set(cell_labels)}
    soft = soft_shapes(poly, cell_labels, stiff_ratio)
    nb = lambda p: math.comb(p + 3, 3)
    n_soft = sum(counts[l] for l in soft)
    for p_big in (6, 5, 4, 3, 2):
        for p_small in (2, 1):
            n = sum(counts[l] * nb(p_big if l in soft else p_small) for l in counts)
            if n <= max_dofs:
                return {l: p_big if l in soft else p_small for l in counts}, n_soft
    return {l: 2 if l in soft else 1 for l in counts}, n_soft


def isolated_walls(poly, cell_labels, stiff_ratio=3.0):
    """(soft shape, one soft cell, its neighbour cells) if there is a single soft shape and
    no two soft walls share a ridge, else None."""
    soft = soft_shapes(poly, cell_labels, stiff_ratio)
    if len(soft) != 1:
        return None
    lab = next(iter(soft))
    owner = {}
    for c, fs in enumerate(poly.cell_faces):
        for f in fs:
            owner.setdefault(f, []).append(c)
    if any(cell_labels[a] == lab and cell_labels[b] == lab for a, b in owner.values()):
        return None
    c0 = cell_labels.index(lab)
    nbrs = sorted({x for f in poly.cell_faces[c0] for x in owner[f] if x != c0})
    return lab, c0, nbrs


def run_isolated(poly, args, cell_class_info, iso, labels):
    """Largest walls that touch only much stiffer walls: every band is one wall mode of all
    soft walls (multiplicity = walls x wall multiplet). Computed on one soft wall welded to all
    its neighbours. With their outer ridges clamped, all other soft walls are at rest: that is
    the band's mean (trace of the band, to first order; checked on prix against the full
    coupled model, 0.4 %). Hinged towards the other soft walls it lies in the lower part of the
    band: the difference estimates the band's half width and is part of the error. Exported up
    to the lowest simply supported fundamental of the neighbour shapes, below which no band of
    the neighbours can lie. Returns None if that leaves fewer than --isolated-min x f1."""
    lab, c0, nbrs = iso
    p, lam_s = args.wall_degree, 2 / 3
    t0 = time.time()
    ss = min(CellShape(poly, labels.index(l), 8, lam_s, 1.0).clamped_omegas(
        poly, labels.index(l), 0.0, H_REF)[0] for l in {labels[c] for c in nbrs})
    single = {}
    for q in (p - 2, p):
        w = CellShape(poly, c0, q, lam_s, 1.0, clamped=True).clamped_omegas(poly, c0, 0.0, H_REF)
        single[q] = [(w[g].mean(), len(g)) for g in group_multiplets(w, 1e-5)]
        if ss < args.isolated_min * single[q][0][0]:
            return None
    nb_labels = {labels[c] for c in nbrs}
    soft_shape = CellShape(poly, c0, p, lam_s, 1.0, p_slope=p - 2)
    patch = {}
    for pn in (args.neighbour_degree - 2, args.neighbour_degree):
        degree = {l: ((p - 2, p) if l == lab else (pn, pn)) for l in nb_labels | {lab}}
        shapes = {l: CellShape(poly, labels.index(l), pn, lam_s, 1.0, p_slope=pn) for l in nb_labels}
        shapes[lab] = soft_shape
        for key, clamp in (("clamped", lambda c: True), ("hinged", lambda c: labels[c] != lab)):
            if pn < args.neighbour_degree and key == "hinged":
                continue
            om, wsq, _, _ = analyse_hollow(poly, labels, degree, args.beta, 300,
                                           cells=[c0] + nbrs, clamp_outside=clamp, shapes=shapes,
                                           verbose=False)
            share = wsq[0] / wsq.sum(axis=0)
            groups = [g for g in group_multiplets(om, 1e-5) if share[g].mean() >= 0.5]
            patch[pn, key] = [(om[g].mean(), len(g)) for g in groups]
            if pn == args.neighbour_degree and key == "clamped":
                spill = [wsq[:, g].sum(axis=1) for g in groups]   # per patch cell, multiplet
    pn = args.neighbour_degree
    n_soft = labels.count(lab)
    freq, mult, err = [], [], []
    for i, (w, m) in enumerate(single[p]):
        if w >= ss or i >= min(len(v) for v in patch.values()):
            break
        (wc, mc), (wh, mh_), (wl, ml) = patch[pn, "clamped"][i], patch[pn, "hinged"][i], \
            patch[pn - 2, "clamped"][i]
        if not mc == mh_ == ml == m:
            raise RuntimeError(f"{poly.name}: patch multiplet {i} does not match the single wall")
        disc = abs(single[p - 2][i][0] - w) / w if i < len(single[p - 2]) else 1.0
        freq.append(float(wc / (2 * np.pi)))
        mult.append(n_soft * m)
        err.append(float(max(disc, abs(wl - wc) / wc, (wc - wh) / wc)))
    # gains: a band sums, per soft wall, its m local modes; each local mode spills a little
    # into the neighbours (from the patch), so a neighbour cell collects the spill of all soft
    # walls it touches. Cells touching no soft wall get 0.
    vol = {l: CellShape(poly, labels.index(l), 1, lam_s, 1.0, bubble="tanh").volume
              * cell_frame(poly, labels.index(l))[2] ** 3 for l in set(labels)}
    area = sum(vol[l] for l in labels)
    mass = H_REF * area
    owner = {}
    for c, fs in enumerate(poly.cell_faces):
        for f in fs:
            owner.setdefault(f, []).append(c)
    patch_cells = [c0] + nbrs
    cell_class, class_labels = cell_class_info
    cells = []
    for s, cl in enumerate(class_labels):
        members = [c for c in range(len(labels)) if cell_class[c] == s]
        v = vol[labels[members[0]]]
        if labels[members[0]] == lab:
            gain = [spill[i][0] / v * mass for i in range(len(mult))]
        else:
            touching = sum(labels[x] == lab for f in poly.cell_faces[members[0]]
                           for x in owner[f] if x != members[0])
            idx = [j for j, c in enumerate(patch_cells) if cell_class[c] == s]
            gain = [touching * np.mean([spill[i][j] for j in idx]) / v * mass if idx else 0.0
                    for i in range(len(mult))]
        cells.append({"label": cl, "count": len(members), "normal": [float(x) for x in gain]})
    print(f"  {poly.name:9s} isolated walls: {n_soft} x {lab}, single wall + {len(nbrs)} neighbours, "
          f"degree {p}; neighbours' lowest simply supported mode {ss / single[p][0][0]:.2f} x f1 "
          f"[{time.time() - t0:.0f}s]", flush=True)
    return {
        "freq": freq, "mult": mult, "err": err, "mean": [float(m) for m in mult],
        "cells": cells, "cell_class": cell_class, "area": float(area),
        "degree": {l: (p if l == lab else pn) for l in nb_labels | {lab}},
        "method": "isolated walls: the largest walls share no ridge and touch only much stiffer "
                  "walls, so each band is one wall mode of all of them (one multiplet per band, "
                  "multiplicity = walls x wall multiplet). Ritz on one such wall welded to all "
                  "its neighbours with their outer ridges clamped (the other large walls at "
                  "rest = the band mean), smooth-distance bubble phi = (sum_f d_f^-2)^(-1/2): "
                  "phi^2 x poly (degree p) + phi x poly (degree p - 2; neighbours p, p). Bands "
                  "up to the lowest simply supported fundamental of the neighbour walls (no "
                  f"neighbour band below it); hinge penalty beta = {args.beta:g} mu h^3 / l",
        "rel_error": "max of: change of the single wall from degree p - 2 to p, change from "
                     "neighbour degree p - 2 to p, and the band's estimated half width (outer "
                     "ridges hinged towards the other large walls instead of clamped)",
    }


def gains_from(poly, labels, shapes, wsq, groups, cell_class_info):
    """Gains per multiplet (mass-relative, like modal_output), normal deflection only: the
    random-hit mean (= multiplicity, every mode bends all its mass) and per cell class."""
    C = len(poly.cells)
    vol = np.array([shapes[labels[c]].r ** 3 * shapes[labels[c]].volume for c in range(C)])
    area = vol.sum()                                   # 3D measure of the boundary
    mass = H_REF * area
    per_mode_mean = wsq.sum(axis=0) / area * mass
    cell_class, class_labels = cell_class_info
    cells = []
    for s, lab in enumerate(class_labels):
        members = [c for c in range(C) if cell_class[c] == s]
        mode_gain = np.mean([wsq[c] / vol[c] for c in members], axis=0) * mass
        cells.append({"label": lab, "count": len(members),
                      "normal": [float(mode_gain[g].sum()) for g in groups]})
    return [float(per_mode_mean[g].sum()) for g in groups], cells, float(area)


def run(poly, args, cell_class_info):
    if args.bubble == "tanh":
        return run_tanh(poly, args, cell_class_info)
    labels = [poly.cell_label(c) for c in range(len(poly.cells))]
    if args.isolated:
        iso = isolated_walls(poly, labels)
        res = run_isolated(poly, args, cell_class_info, iso, labels) if iso else None
        if res is not None:
            return res
    return run_dist(poly, args, cell_class_info, labels)


def basis_size(spec):
    """Number of basis functions of a CellShape spec (bubble, slope degree, interior degree)."""
    kind, ps, p = spec
    nb = lambda q: math.comb(q + 3, 3)
    return (2 if kind == "both" else 1) * nb(ps) + nb(p)


def run_dist(poly, args, cell_class_info, labels):
    """Coupled model, block-diagonalised by symmetry. Soft walls: both slope bubbles, slopes
    p - 2, interior p. Stiff walls whose simply supported fundamental lies in the computed
    range take part in the modes: both slope bubbles, slopes --stiff-degree, interior 4; the
    others only transmit rotations: phi slopes --stiff-degree, interior 2. If a symmetry
    block would exceed --block-budget unknowns, the range is lowered (by 0.5 x f1, down to
    3 x f1). All modes up to that range x f1 (at most --budget), exported up to the end of
    the band that reaches over it. Error indicator: change of two symmetry blocks when every
    degree is lowered by 2 (pessimistic)."""
    soft = soft_shapes(poly, labels)
    ps, pst = args.soft_degree, args.stiff_degree
    lam_s = 2 / 3
    counts = {l: labels.count(l) for l in set(labels)}
    rep = {l: labels.index(l) for l in counts}
    made = {}

    def make(l, spec):
        if (l, spec) not in made:
            made[l, spec] = CellShape(poly, rep[l], spec[2], lam_s, 1.0, p_slope=spec[1],
                                      bubble=spec[0])
        return made[l, spec]

    soft_spec = ("both", ps - 2, ps)
    big = max(soft, key=lambda l: make(l, soft_spec).r * min(make(l, soft_spec).d))
    w1 = make(big, soft_spec).clamped_omegas(poly, rep[big], args.beta, H_REF)[0]
    quasi = ("dist", pst, 2)
    w_ss = {l: make(l, quasi).clamped_omegas(poly, rep[l], 0.0, H_REF)[0]
            for l in counts if l not in soft}
    nG = len(mirror_group(poly)[1])
    ratio = args.max_ratio
    while True:
        # upper bound of f1: the largest wall clamped; computed 20 % beyond the export
        # range, so that a band reaching over ratio x f1 can be completed
        omega_max = 1.2 * ratio * w1
        spec = {l: soft_spec if l in soft else
                ("both", pst, 4) if w_ss[l] < omega_max else quasi for l in counts}
        per_block = sum(counts[l] * basis_size(spec[l]) for l in counts) / nG
        if per_block <= args.block_budget or ratio <= 3.0:
            break
        ratio -= 0.5
    shapes = {l: make(l, spec[l]) for l in counts}
    # expected modes below omega_max: per wall, between its simply supported and clamped count
    n_exp = 0.0
    for l, sh in shapes.items():
        w_ss = np.sqrt(np.maximum(np.linalg.eigvalsh(sh.kb), 0) * H_REF ** 2 / 12 / sh.r ** 4)
        w_cl = sh.clamped_omegas(poly, rep[l], args.beta, H_REF)
        n_exp += counts[l] * math.sqrt((w_ss < omega_max).sum() * max((w_cl < omega_max).sum(), 1))
    density = n_exp / sum(counts[l] * shapes[l].size for l in counts)
    omega, wsq, complete, shapes, chi, nG = analyse_hollow_sym(
        poly, labels, spec, args.beta, omega_max, args.budget, shapes=shapes, density=density)
    # export: all multiplets up to max-ratio x f1, plus the rest of a band reaching over it
    # (up to the next gap >= 3 %); if that band is not complete in the computed range (budget),
    # end at the last gap >= 3 % below instead
    gs = group_multiplets(omega, 1e-5)
    fr = np.array([omega[g].mean() for g in gs])
    lim = min(complete, omega_max)
    top = ratio * fr[0]
    if lim >= 1.1 * top:
        # all computed up to 1.1 x top: every multiplet up to top, then the rest of a band
        # reaching over it, up to the next gap >= 3 % (a dense spectrum: at most 1.1 x top)
        n = int(np.searchsorted(fr, top))
        while 0 < n < len(gs) and fr[n] <= 1.1 * top and (fr[n] - fr[n - 1]) / fr[n - 1] < 0.03:
            n += 1
    else:
        # the computed range ends early (mode budget): end at the last gap >= 3 % below it
        n = int(np.searchsorted(fr, lim * (1 - 1e-9)))
        big_gaps = np.flatnonzero((fr[1:n] - fr[:n - 1]) / fr[:n - 1] >= 0.03)
        n = int(big_gaps[-1]) + 1 if len(big_gaps) else n
    groups = gs[:n]
    # error indicator: the same two blocks with every degree lowered by 2
    coarse = {l: (k, a - 2, b - 2 if b > 2 else 2) for l, (k, a, b) in spec.items()}
    err_chis = sorted({0, nG - 1})[:args.err_blocks]
    of, _, _, _, chif, _ = analyse_hollow_sym(
        poly, labels, coarse, args.beta, 1.3 * omega_max, args.budget, characters=err_chis,
        shapes={l: make(l, coarse[l]) for l in counts}, verbose=False)
    rel = np.full(len(omega), np.nan)
    for s in err_chis:
        idx = np.flatnonzero(chi == s)
        fw = of[chif == s]
        m = min(len(idx), len(fw))
        rel[idx[:m]] = np.abs(fw[:m] - omega[idx[:m]]) / omega[idx[:m]]
    known = np.flatnonzero(~np.isnan(rel))
    err = []
    for g in groups:
        e = rel[g][~np.isnan(rel[g])]
        if len(e) == 0 and len(known):                  # no member in those blocks: nearest one
            e = rel[known[np.argmin(np.abs(omega[known] - omega[g[0]]))]]
        err.append(float(np.max(e)) if np.size(e) else float("nan"))
    mean, cells, area = gains_from(poly, labels, shapes, wsq, groups, cell_class_info)
    return {
        "freq": [float(omega[g].mean() / (2 * np.pi)) for g in groups],
        "mult": [len(g) for g in groups], "err": err, "mean": mean, "cells": cells,
        "cell_class": cell_class_info[0], "area": area,
        "degree": {l: (s[2] if l in soft else s[1]) for l, s in spec.items()},
        "basis": {l: f"{k}: slopes {a}, interior {b}" for l, (k, a, b) in sorted(spec.items())},
        "range": ratio,
        "method": "Ritz per cell with two bubbles: phi = (sum_f d_f^-2)^(-1/2), a smooth distance "
                  "to the face planes, and b = prod_f tanh(d_f/delta); phi^2 x poly (zero slope "
                  "at the faces: the interior) plus phi x poly and b x poly for the hinge "
                  "rotations (near an edge they behave like r and r^2, the exact solution like "
                  "r^(pi/alpha) in between). Soft walls: slopes p - 2, interior p; stiff walls "
                  f"resonating in the range: slopes {pst}, interior 4; the others: phi x poly "
                  f"slopes {pst}, interior 2 (see basis; polynomial_degree lists p for soft "
                  "walls, the slope degree for stiff ones). Legendre-product polynomials, "
                  "orthonormalised by QR; cell matrices per shape mapped by congruence; problem "
                  f"block-diagonalised by (Z2)^{int(math.log2(nG))} mirror symmetries (one "
                  "block per character, assembled from one representative cell per orbit); "
                  f"hinge penalty beta = {args.beta:g} mu h^3 / l",
        "rel_error": "indicator: change of the modes of two symmetry blocks when every degree "
                     "is lowered by 2 (per multiplet the largest change of its members, else "
                     "of the nearest mode); pessimistic: the step from these degrees to 2 "
                     "higher is several times smaller. Ritz values are upper bounds",
    }


def run_tanh(poly, args, cell_class_info):
    """First version: flat-top tanh bubble, no symmetry reduction (reproduces the files of
    before 2026-09-27; too high by up to ~30 % for walls with many faces)."""
    labels = [poly.cell_label(c) for c in range(len(poly.cells))]
    iso = isolated_walls(poly, labels) if args.isolated else None
    if iso is not None:
        res = run_isolated(poly, args, cell_class_info, iso, labels)
        if res is not None:
            return res
    degree, n_soft = choose_degrees(poly, labels, args.max_dofs)
    p_big, p_small = max(degree.values()), min(degree.values())
    # enough modes for the first bands of the soft walls (a band holds one mode per soft wall
    # and wall mode; isolated equal walls give very narrow, highly degenerate bands)
    n_modes = int(min(max(args.modes, 6 * n_soft), args.max_modes))
    omega, wsq, degree, shapes = analyse_hollow(poly, labels, degree, args.beta, n_modes,
                                                bubble="tanh", delta_frac=args.delta)
    # error indicator: per soft cell shape, the relative change of its clamped-wall
    # frequencies (single wall, welds' slope penalty on all faces) from degree p to p + 2,
    # taken at the wall mode nearest to each global mode and weighted with the mode's share
    # of bending energy in that shape
    lam_s = 2 / 3
    share = {lab: np.zeros(len(omega)) for lab in degree}
    for c, lab in enumerate(labels):
        share[lab] += wsq[c]
    total = sum(share.values())
    rel_err = np.zeros(len(omega))
    for lab, p in degree.items():
        if p == p_small and p != p_big:
            continue                                   # stiff walls only transmit rotations
        sh = shapes[lab]
        fine = CellShape(poly, labels.index(lab), p + 2, lam_s, 1.0, delta_frac=args.delta,
                         bubble="tanh")
        rep = labels.index(lab)
        wa = sh.clamped_omegas(poly, rep, args.beta, H_REF)
        wb = fine.clamped_omegas(poly, rep, args.beta, H_REF)
        m = min(len(wa), len(wb))
        diff = np.abs(wa[:m] - wb[:m]) / wb[:m]
        nearest = np.abs(wa[None, :m] - omega[:, None]).argmin(axis=1)
        rel_err += share[lab] / total * diff[nearest]
    # export range: up to --max-ratio x the lowest frequency, ending with a complete band
    groups = []
    for g in group_multiplets(omega, 1e-5)[:-1]:      # last one may be cut
        if omega[g[0]] > args.max_ratio * omega[0]:
            break
        groups.append(g)
    groups = groups[:complete_bands([omega[g[0]] for g in groups])]
    mean, cells, area = gains_from(poly, labels, shapes, wsq, groups, cell_class_info)
    return {
        "freq": [float(omega[g].mean() / (2 * np.pi)) for g in groups],
        "mult": [len(g) for g in groups],
        "err": [float(rel_err[g].max()) for g in groups],
        "mean": mean, "cells": cells, "cell_class": cell_class_info[0], "degree": degree,
        "area": area,
    }


def write_output(poly, res, out_dir, args):
    zeros = [0.0] * len(res["freq"])
    gains = lambda normal: {"normal": [sig(x, 5) for x in normal], "tangential": zeros}
    doc = {
        "name": poly.name,
        "description": poly.description,
        "source": os.path.relpath(poly.path).replace("\\", "/"),
        "format": FORMAT,
        "model": {
            "body": "thin-walled hollow 4D body: the cells are welded 3D plates of thickness h",
            "equations": "Kirchhoff bending of each cell along its 4D normal, (h^3/12)[lam* (Lap w)^2 "
                         "+ 2 mu |grad grad w|^2], lam* = 2 lam mu/(lam+2mu); ridges fixed "
                         "(thin-wall limit), hinge rotation continuous across ridges (welded)",
            "method": res.get("method",
                              "Ritz per cell: flat-top bubble prod tanh(dist_f/delta) times "
                              "polynomials, cell matrices per shape mapped by congruence, hinge "
                              f"penalty beta = {args.beta:g} mu h^3 / l"),
            "polynomial_degree": {k: v for k, v in sorted(res["degree"].items())},
            **({"basis": res["basis"]} if "basis" in res else {}),
            **({"computed_up_to": f"{res['range']:g} x f1 (plus the band reaching over it)"}
               if "range" in res else {}),
            "lambda_over_mu": 1.0,
            "units": "circumradius R = 1, density rho = 1, shear modulus mu = 1 (shear wave "
                     f"speed c_s = 1); wall thickness h_ref = {H_REF} R",
            "frequency": "f_Hz = frequency * (h / h_ref) * c_s[m/s] / R[m]  (bending: f ~ h)",
            "gains": "as in modal_output (mass-relative, per multiplet sum); mean = multiplicity "
                     "(every mode bends all its mass); tangential = 0 and no vertex gains "
                     "(membrane motion is not modelled, ridges and vertices are fixed)",
            "rel_error": res.get("rel_error",
                                 "indicator: change of the clamped single-wall frequencies of the "
                                 "dominant cell shapes from degree p to p + 2; convergence is slow "
                                 "for cells with obtuse inner dihedral angles, frequencies are "
                                 "upper bounds"),
            "multiplet_merge_tol": 1e-5,
        },
        "geometry": {
            "boundary_volume": sig(res["area"], 10),
            "edge_length": sig(1 / poly.radius, 10),
            "counts": [len(poly.xi), len(poly.edges), len(poly.faces), len(poly.cells)],
        },
        "modes": {
            "frequency": [sig(f, 8) for f in res["freq"]],
            "multiplicity": res["mult"],
            "rel_error": [sig(e, 2) for e in res["err"]],
        },
        "gains": {
            "mean": gains(res["mean"]),
            "cells": [dict(label=s["label"], count=s["count"], **gains(s["normal"]))
                      for s in res["cells"]],
        },
        "cell_class": res["cell_class"],
    }
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, poly.name + ".json")
    with open(path + ".tmp", "w", encoding="utf-8") as f:
        json.dump(doc, f, separators=(",", ":"))
    os.replace(path + ".tmp", path)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("names", nargs="*", help="polychora to analyse (default: all)")
    ap.add_argument("--out", default="modal_hollow_output")
    ap.add_argument("--solid", default="modal_output",
                    help="solid-body results: their cell classes (symmetry orbits) are reused")
    ap.add_argument("--max-ratio", type=float, default=5.0,
                    help="export multiplets up to this multiple of the lowest frequency")
    ap.add_argument("--budget", type=int, default=20000,
                    help="at most about this many modes per polychoron (fewer: lower range)")
    ap.add_argument("--soft-degree", type=int, default=8,
                    help="interior degree of the soft (large) walls; their slopes: 2 less")
    ap.add_argument("--stiff-degree", type=int, default=6,
                    help="slope degree of the stiff (small) walls; interior 4 if they resonate in the "
                         "range, else 2")
    ap.add_argument("--block-budget", type=int, default=32000,
                    help="largest symmetry block (unknowns); above it the range is lowered")
    ap.add_argument("--err-blocks", type=int, default=2,
                    help="symmetry blocks recomputed with degrees - 2 for the error indicator")
    ap.add_argument("--beta", type=float, default=1e3, help="hinge penalty factor")
    ap.add_argument("--isolated", action="store_true",
                    help="isolated largest walls (prahi, prix, gidpixhi): one wall + neighbours "
                         "(band means only) instead of the coupled model")
    ap.add_argument("--isolated-min", type=float, default=3.5,
                    help="isolated-wall model only if the neighbours' lowest simply supported "
                         "mode is at least this multiple of f1 (else the coupled model)")
    ap.add_argument("--wall-degree", type=int, default=10,
                    help="polynomial degree of the isolated-wall model ...")
    ap.add_argument("--neighbour-degree", type=int, default=8, help="... and of its neighbours")
    ap.add_argument("--bubble", choices=["dist", "tanh"], default="dist",
                    help="dist: the current model (README); tanh: the first version (reproduce it "
                         "with --bubble tanh --max-ratio 3)")
    ap.add_argument("--modes", type=int, default=400,
                    help="tanh: modes to compute, at least 6 per wall of the most frequent shape ...")
    ap.add_argument("--max-modes", type=int, default=800, help="tanh: ... up to this many")
    ap.add_argument("--max-dofs", type=int, default=26000,
                    help="tanh: budget that sets the polynomial degrees per polychoron")
    ap.add_argument("--delta", type=float, default=None,
                    help="tanh: bubble layer delta / inradius (default: 3 for cells with <= 8 "
                         "faces, 1 for <= 12, 0.5 above)")
    args = ap.parse_args()

    paths = sorted(glob.glob(os.path.join("topology_output", "*.json")))
    if args.names:
        paths = [os.path.join("topology_output", n + ".json") for n in args.names]
    print(f"hollow modal analysis: h_ref = {H_REF} R, {args.bubble} basis, up to "
          f"{args.max_ratio:g} f1, {len(paths)} polychora -> {args.out}/", flush=True)
    def n_cells(p):
        with open(p, encoding="utf-8-sig") as f:
            return len(json.load(f)["cells"])
    paths.sort(key=n_cells)                             # small ones first
    failed = []
    for p in paths:
        t0 = time.time()
        try:
            poly = Polytope(p)
            with open(os.path.join(args.solid, poly.name + ".json"), encoding="utf-8") as f:
                solid = json.load(f)
            info = (solid["cell_class"], [c["label"] for c in solid["gains"]["cells"]])
            res = run(poly, args, info)
            write_output(poly, res, args.out, args)
        except Exception as e:                          # keep going with the others
            import traceback
            traceback.print_exc()
            failed.append(os.path.basename(p))
            print(f"  FAILED {p}: {e}", flush=True)
            continue
        f = np.array(res["freq"])
        print(f"      {sum(res['mult'])} modes / {len(f)} multiplets up to {f[-1] / f[0]:.2f} f1, "
              f"max rel. error {max(res['err']):.1%} (first multiplet {res['err'][0]:.1%})  "
              + " ".join(f"{x / f[0]:.3f}x{m}" for x, m in list(zip(f, res["mult"]))[:6])
              + f"  [{(time.time() - t0) / 60:.1f} min]", flush=True)
    if failed:
        print("failed:", " ".join(failed))
        sys.exit(1)


if __name__ == "__main__":
    main()

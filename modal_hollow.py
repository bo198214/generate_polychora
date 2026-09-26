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

Discretisation
--------------
Per cell, w = b(s) * poly(s) with poly of total degree <= p in local 3D coordinates and a
"flat-top" bubble b = prod_faces tanh(dist_f / delta) (delta = inradius / 3): it vanishes
linearly on every face like the polynomial bubble prod dist_f (its delta -> infinity limit),
but stays ~1 inside even for cells with 62 faces. The limit p -> infinity does not depend on
delta. Cell matrices are computed once per cell shape (flag-tetrahedra Gauss quadrature) and
carried to every congruent cell by an orthogonal map found from the vertex sets; cells whose
walls are much stiffer than the softest ones get a lower degree (they only transmit rotations
in the audible band). Error indicator per mode: change of the dominant cell shapes'
clamped single-wall frequencies from degree p to p + 2.

Units: circumradius R = 1, rho = 1, mu = 1, lambda = mu (as modal_output), reference wall
thickness h_ref = 0.06 R. f_Hz = frequency * (h / h_ref) * c_s / R.
"""
import argparse
import glob
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
    """Reference cell of one shape: face planes, flat-top bubble basis, bending/mass matrices.
    Local coordinates s = T (x - o) / r (r = cell circumradius); matrices are per unit
    thickness factors: K = kb * h^3 / 12 / r, M = G * h * r^3 with the returned kb, G."""

    def __init__(self, poly, c, p, lam_s, mu, q=None, delta_frac=1 / 3):
        self.o, self.T, self.r = cell_frame(poly, c)
        self.p = p
        self.mono = Monomials(3, p)
        self.E = self.mono.upto(p)
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
        self.delta = self.d.min() * delta_frac
        self.norm = np.prod(np.tanh(self.d / self.delta))
        E = self.E
        pos = {tuple(e): k for k, e in enumerate(E)}
        self.dmap = []
        for i in range(3):
            ei = np.eye(3, dtype=int)[i]
            self.dmap.append((np.array([pos.get(tuple(e - ei), 0) for e in E]), E[:, i].astype(float)))
        # quadrature on the flag tetrahedra (centroid, face centre, edge endpoints)
        qq = q or p + 6
        Z, W = tet_rule(qq)
        S, Wt = [], []
        for f in poly.cell_faces[c]:
            for e in poly.face_edges[f]:
                a, bv = poly.xi[poly.edges[e]]
                Vt = loc(np.stack([self.o, poly.face_ctr[f], a, bv]))
                D = Vt[1:] - Vt[:-1]
                S.append(Vt[0] + Z @ D)
                Wt.append(W * abs(np.linalg.det(D)))
        S, Wt = np.concatenate(S), np.concatenate(Wt)
        nb = len(E)
        kb, G = np.zeros((nb, nb)), np.zeros((nb, nb))
        for s0 in range(0, len(S), 4000):
            Ws = Wt[s0:s0 + 4000]
            phi, _, hess = self.eval(S[s0:s0 + 4000])
            lap = hess[0, 0] + hess[1, 1] + hess[2, 2]
            kb += lam_s * lap.T @ (Ws[:, None] * lap)
            for i in range(3):
                for j in range(3):
                    kb += 2 * mu * hess[i, j].T @ (Ws[:, None] * hess[i, j])
            G += phi.T @ (Ws[:, None] * phi)
        self.kb, self.G, self.volume = (kb + kb.T) / 2, (G + G.T) / 2, float(Wt.sum())

    def local(self, X):
        return (X - self.o) @ self.T.T / self.r

    def _mono(self, s):
        mono = self.mono.evaluate(s, self.p)[:, :len(self.E)]
        dmono = [mono[:, self.dmap[i][0]] * self.dmap[i][1] for i in range(3)]
        return mono, dmono

    def eval(self, s):
        """phi (m, nb), grad (3, m, nb), hessian (3, 3, m, nb) at interior points s (local)."""
        L = self.d[None, :] - s @ self.m.T
        t = np.tanh(L / self.delta)
        sech2 = 1 - t * t
        g = sech2 / (self.delta * t)                    # t'/t
        c2 = -2 * sech2 / self.delta ** 2 - g * g       # t''/t - (t'/t)^2
        b = np.prod(t, axis=1) / self.norm
        gv = g @ self.m                                  # sum_f g_f m_f   (grad b = -b gv)
        gb = -b[:, None] * gv
        Hb = b[None, None, :] * (np.einsum("mi,mj->ijm", gv, gv)
                                 + np.einsum("mf,fi,fj->ijm", c2, self.m, self.m))
        mono, dmono = self._mono(s)
        phi = b[:, None] * mono
        grad = np.stack([gb[:, i:i + 1] * mono + b[:, None] * dmono[i] for i in range(3)])
        hess = np.empty((3, 3, len(s), len(self.E)))
        for i in range(3):
            for j in range(3):
                d2 = dmono[i][:, self.dmap[j][0]] * self.dmap[j][1]
                hess[i, j] = (Hb[i, j][:, None] * mono + gb[:, i:i + 1] * dmono[j]
                              + gb[:, j:j + 1] * dmono[i] + b[:, None] * d2)
        return phi, grad, hess

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
        L = self.d[None, :] - s @ self.m.T
        L[:, k] = 0.0
        t = np.tanh(L / self.delta)
        t[:, k] = 1.0
        db = np.prod(t, axis=1) / (self.delta * self.norm)
        mono, _ = self._mono(s)
        return db[:, None] * mono


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


def analyse_hollow(poly, cell_labels, p_big, p_small, beta, n_modes, lam=1.0, mu=1.0,
                   stiff_ratio=3.0, delta_frac=1 / 3, verbose=True):
    """Welded thin-walled hollow polychoron. Returns frequencies (h = H_REF), multiplets,
    per-cell ∫ w^2 contributions for gains."""
    t0 = time.time()
    lam_s = 2 * lam * mu / (lam + 2 * mu)
    h = H_REF
    C = len(poly.cells)
    labels = sorted(set(cell_labels))
    # which shapes are soft (in the audible band) -> degree p_big, else p_small
    inr = {}
    for lab in labels:
        c = cell_labels.index(lab)
        o, T, r = cell_frame(poly, c)
        fdist = [np.linalg.norm(poly.face_ctr[f] - o) for f in poly.cell_faces[c]]
        inr[lab] = min(fdist)
    soft = max(inr.values())
    degree = {lab: (p_big if (soft / inr[lab]) ** 2 < stiff_ratio else p_small) for lab in labels}
    shapes = {}
    for lab in labels:
        rep = cell_labels.index(lab)
        shapes[lab] = CellShape(poly, rep, degree[lab], lam_s, mu, delta_frac=delta_frac)
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
    Qs, frames_ = [], []
    ref_adj = {lab: cell_adj(cell_labels.index(lab)) for lab in labels}
    for c in range(C):
        sh = shapes[cell_labels[c]]
        o, T, r = cell_frame(poly, c)
        V = (poly.xi[poly.cells[c]] - o) @ T.T / r
        Q = congruence(sh, V, ref_adj[cell_labels[c]], cell_adj(c))
        if Q is None:
            raise RuntimeError(f"{poly.name}: cell {c} not congruent to its reference")
        Qs.append(Q)
        frames_.append((o, T, r))
    t1 = time.time()

    sizes = [len(shapes[cell_labels[c]].E) for c in range(C)]
    offs = np.concatenate([[0], np.cumsum(sizes)])
    N = int(offs[-1])
    Kb, Mb = {}, {}
    for c in range(C):
        sh = shapes[cell_labels[c]]
        Kb[c, c] = h ** 3 / 12 / sh.r * sh.kb
        Mb[c, c] = h * sh.r ** 3 * sh.G
    ell = float(np.mean([shapes[l].r for l in labels]))
    kr = beta * mu * h ** 3 / ell

    face_cells = [[] for _ in poly.faces]
    for c, fs in enumerate(poly.cell_faces):
        for f in fs:
            face_cells[f].append(c)
    U, W = tri_rule(max(degree.values()) + 5)
    for f, (A, B) in enumerate(face_cells):
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
        for Xc in (A, B):
            sh = shapes[cell_labels[Xc]]
            o, T, r = frames_[Xc]
            s = ((X - o) @ T.T / r) @ Qs[Xc]           # reference-cell coordinates
            k = list(poly.cell_faces[Xc]).index(f)
            # the same face in the reference cell: the reference face whose plane contains s
            mref = s @ sh.m.T - sh.d[None, :]
            kf = int(np.argmin(np.abs(mref).max(axis=0)))
            gs.append(sh.face_slope(kf, s) / r)
            m_in = -(Qs[Xc] @ sh.m[kf])                  # inward, in the cell's own local frame
            nm.append((poly.cell_normal[Xc], T.T @ m_in))
        (nA, mA), (nB, mB) = nm
        sgn = np.sign(nA @ nB * (mA @ mB) - nA @ mB * (mA @ nB))
        gA, gB = gs[0], -sgn * gs[1]
        for I, gI in ((A, gA), (B, gB)):
            for J, gJ in ((A, gA), (B, gB)):
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
    wsq = np.zeros((C, len(omega)))
    for c in range(C):
        Xc = X[offs[c]:offs[c + 1]]
        wsq[c] = np.einsum("ak,ak->k", Xc, Mb[c, c] @ Xc) / h
    if verbose:
        print(f"  {poly.name:9s} cells={C:5d} DOFs={N:6d} degrees="
              + ",".join(f"{lab.split()[-1][:4]}:{degree[lab]}" for lab in labels)
              + f"  f1={omega[0] / (2 * np.pi):.5f}  "
              f"[{t1 - t0:.0f}s shapes, {t2 - t1:.0f}s coupling, {t3 - t2:.0f}s eigen]", flush=True)
    return omega, wsq, degree, shapes


def choose_degrees(poly, cell_labels, max_dofs, stiff_ratio=3.0):
    """Highest degree for the soft cells (and 2 or 1 for stiff ones) within the DOF budget;
    also returns how many walls are soft."""
    counts = {lab: cell_labels.count(lab) for lab in set(cell_labels)}
    inr = {}
    for lab in counts:
        c = cell_labels.index(lab)
        o = poly.xi[poly.cells[c]].mean(axis=0)
        inr[lab] = min(np.linalg.norm(poly.face_ctr[f] - o) for f in poly.cell_faces[c])
    soft = max(inr.values())
    is_soft = {lab: (soft / inr[lab]) ** 2 < stiff_ratio for lab in counts}
    nb = lambda p: math.comb(p + 3, 3)
    for p_big in (6, 5, 4, 3, 2):
        for p_small in (2, 1):
            n = sum(counts[l] * nb(p_big if is_soft[l] else p_small) for l in counts)
            if n <= max_dofs:
                return p_big, p_small, sum(counts[l] for l in counts if is_soft[l])
    return 2, 1, sum(counts[l] for l in counts if is_soft[l])


def run(poly, args, cell_class_info):
    labels = [poly.cell_label(c) for c in range(len(poly.cells))]
    p_big, p_small, n_soft = choose_degrees(poly, labels, args.max_dofs)
    # enough modes for the first bands of the soft walls (a band holds one mode per soft wall
    # and wall mode; isolated equal walls give very narrow, highly degenerate bands)
    n_modes = int(min(max(args.modes, 6 * n_soft), args.max_modes))
    omega, wsq, degree, shapes = analyse_hollow(poly, labels, p_big, p_small, args.beta, n_modes,
                                                delta_frac=args.delta)
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
        fine = CellShape(poly, labels.index(lab), p + 2, lam_s, 1.0, delta_frac=args.delta)
        rep = labels.index(lab)
        wa = sh.clamped_omegas(poly, rep, args.beta, H_REF)
        wb = fine.clamped_omegas(poly, rep, args.beta, H_REF)
        m = min(len(wa), len(wb))
        diff = np.abs(wa[:m] - wb[:m]) / wb[:m]
        nearest = np.abs(wa[None, :m] - omega[:, None]).argmin(axis=1)
        rel_err += share[lab] / total * diff[nearest]
    # export range: up to --max-ratio x the lowest frequency
    groups = []
    for g in group_multiplets(omega, 1e-5)[:-1]:      # last one may be cut
        if omega[g[0]] > args.max_ratio * omega[0]:
            break
        groups.append(g)

    # gains (mass-relative, like modal_output): normal deflection only (bending model)
    C = len(poly.cells)
    vol = np.array([shapes[labels[c]].r ** 3 * shapes[labels[c]].volume for c in range(C)])
    area = vol.sum()                                   # 3D measure of the boundary
    mass = H_REF * area
    per_mode_mean = wsq.sum(axis=0) / area * mass      # = 1 for every mode (all mass bends)
    cell_class, class_labels = cell_class_info
    cells = []
    for s, lab in enumerate(class_labels):
        members = [c for c in range(C) if cell_class[c] == s]
        mode_gain = np.mean([wsq[c] / vol[c] for c in members], axis=0) * mass
        cells.append({"label": lab, "count": len(members),
                      "normal": [float(mode_gain[g].sum()) for g in groups]})
    freq = [float(omega[g].mean() / (2 * np.pi)) for g in groups]
    return {
        "freq": freq, "mult": [len(g) for g in groups],
        "err": [float(rel_err[g].max()) for g in groups],
        "mean": [float(per_mode_mean[g].sum()) for g in groups],
        "cells": cells, "cell_class": cell_class, "degree": degree,
        "area": float(area),
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
            "method": "Ritz per cell: flat-top bubble prod tanh(dist_f/delta) times polynomials, "
                      "cell matrices per shape mapped by congruence, hinge penalty "
                      f"beta = {args.beta:g} mu h^3 / l",
            "polynomial_degree": {k: v for k, v in sorted(res["degree"].items())},
            "lambda_over_mu": 1.0,
            "units": "circumradius R = 1, density rho = 1, shear modulus mu = 1 (shear wave "
                     f"speed c_s = 1); wall thickness h_ref = {H_REF} R",
            "frequency": "f_Hz = frequency * (h / h_ref) * c_s[m/s] / R[m]  (bending: f ~ h)",
            "gains": "as in modal_output (mass-relative, per multiplet sum); mean = multiplicity "
                     "(every mode bends all its mass); tangential = 0 and no vertex gains "
                     "(membrane motion is not modelled, ridges and vertices are fixed)",
            "rel_error": "indicator: change of the clamped single-wall frequencies of the dominant "
                         "cell shapes from degree p to p + 2; convergence is slow for cells with obtuse "
                         "inner dihedral angles, frequencies are upper bounds",
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
    ap.add_argument("--modes", type=int, default=400,
                    help="modes (with multiplicity) to compute, at least 6 per wall of the "
                         "most frequent shape ...")
    ap.add_argument("--max-modes", type=int, default=800, help="... up to this many")
    ap.add_argument("--max-dofs", type=int, default=14000,
                    help="budget that sets the polynomial degrees per polychoron")
    ap.add_argument("--max-ratio", type=float, default=3.0,
                    help="export multiplets up to this multiple of the lowest frequency")
    ap.add_argument("--beta", type=float, default=1e3, help="hinge penalty factor")
    ap.add_argument("--delta", type=float, default=0.5, help="bubble layer delta / inradius")
    args = ap.parse_args()

    paths = sorted(glob.glob(os.path.join("topology_output", "*.json")))
    if args.names:
        paths = [os.path.join("topology_output", n + ".json") for n in args.names]
    print(f"hollow modal analysis: h_ref = {H_REF} R, <= {args.max_dofs} DOFs, "
          f"{len(paths)} polychora -> {args.out}/", flush=True)
    for p in paths:
        poly = Polytope(p)
        with open(os.path.join(args.solid, poly.name + ".json"), encoding="utf-8") as f:
            solid = json.load(f)
        info = (solid["cell_class"], [c["label"] for c in solid["gains"]["cells"]])
        res = run(poly, args, info)
        write_output(poly, res, args.out, args)
        f = np.array(res["freq"])
        print(f"      {sum(res['mult'])} modes / {len(f)} multiplets up to {f[-1] / f[0]:.2f} f1, "
              f"max rel. error {max(res['err']):.1%} (first multiplet {res['err'][0]:.1%})  "
              + " ".join(f"{x / f[0]:.3f}x{m}" for x, m in list(zip(f, res["mult"]))[:6]),
              flush=True)


if __name__ == "__main__":
    main()

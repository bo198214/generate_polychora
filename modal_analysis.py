#!/usr/bin/env python3
"""
Modal analysis of the polychora as vibrating solid 4D bodies.

Usage:
  python modal_analysis.py                  # every polychoron in topology_output/
  python modal_analysis.py pen tes ex       # selected polychora (by name)
  python modal_analysis.py --selftest       # verify against analytic 3D/4D ball solutions

Output: modal_output/<name>.json (format: see README.md, "Modal analysis").

Physics
-------
Each polychoron is a homogeneous, isotropic, linear-elastic *solid* 4D body, free
(unsupported) in R^4. Small vibrations u(x,t) in R^4 obey the 4D Navier-Cauchy equations

    rho u_tt = (lambda + mu) grad(div u) + mu Laplace(u)        inside P
    sigma(u) n = 0                                              on the boundary (3-manifold)

with sigma = lambda tr(eps) I + 2 mu eps, eps = sym(grad u). Normal modes
u = phi(x) sin(omega t) solve  K phi = omega^2 M  phi.  A free 4D body has 10
rigid-body modes (4 translations + 6 rotation planes) at omega = 0; they are removed.

Units: circumradius R = 1, density rho = 1, shear modulus mu = 1, i.e. the shear wave
speed c_s = 1. A dimensionless frequency f therefore means  f_Hz = f * c_s / R.

Discretisation
--------------
Rayleigh-Ritz with a complete polynomial basis: every displacement component is a
polynomial of total degree <= N in (x, y, z, w) -- the 4D version of Visscher's "xyz
algorithm" from resonant ultrasound spectroscopy. Mass and stiffness entries are
integrals of monomials over the polytope; they are computed exactly (up to rounding)
with Lasserre's divergence-theorem recursion vertex -> edge -> face -> cell -> 4-polytope:
for a k-face G, a point x0 in its affine hull and homogeneous f of degree q,

    (k + q) * int_G f = sum_{facets H of G} dist(x0, H) * int_H f  +  int_G x0 . grad f .

Accuracy is estimated by repeating the solve with degree N-2 (Ritz values converge from
above); modes are exported up to the first multiplet whose estimate exceeds --max-error.

Symmetry
--------
Coordinate reflections / central inversion of the polytope are detected and used to
split the Ritz problem into independent parity blocks (exact, just faster). Degenerate
multiplets (from the polytope's symmetry group) are grouped. Mode shapes are exported
only as excitation gains summed over each multiplet (independent of the arbitrary basis
inside a degenerate eigenspace): averaged over each class of symmetry-equivalent cells
(exactly, via per-cell monomial moments), over the whole boundary, and at a vertex.
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
from scipy.spatial import cKDTree

FORMAT = "polychoron-modal/1"


def n_rigid(n):
    """Rigid-body modes of a free body in R^n: n translations + n(n-1)/2 rotation planes."""
    return n + n * (n - 1) // 2


# ─── monomials ────────────────────────────────────────────────────────────

def _compositions(d, n):
    if n == 1:
        yield (d,)
        return
    for k in range(d, -1, -1):
        for rest in _compositions(d - k, n - 1):
            yield (k,) + rest


class Monomials:
    """All monomials x^g in n variables with |g| <= D, grouped by total degree."""

    def __init__(self, n, D):
        self.n, self.D = n, D
        self.levels = [np.array(list(_compositions(d, n)), dtype=np.intp).reshape(-1, n)
                       for d in range(D + 1)]
        self.lut = np.zeros((D + 1,) * n, dtype=np.intp)      # index within own level
        for lev in self.levels:
            self.lut[tuple(lev.T)] = np.arange(len(lev))
        # evaluation: x^g = x_l * x^(g - e_l), l = first axis with g_l > 0
        self.mul_axis, self.mul_pred = [None], [None]
        # Lasserre gradient term: x0 . grad(x^g) = sum_l x0_l g_l x^(g - e_l)
        self.grad = [None]
        for d in range(1, D + 1):
            lev = self.levels[d]
            ax = np.argmax(lev > 0, axis=1)
            prev = lev.copy()
            prev[np.arange(len(lev)), ax] -= 1
            self.mul_axis.append(ax)
            self.mul_pred.append(self.lut[tuple(prev.T)])
            terms = []
            for l in range(n):
                cols = np.flatnonzero(lev[:, l] > 0)
                prev = lev[cols].copy()
                prev[:, l] -= 1
                terms.append((l, cols, self.lut[tuple(prev.T)], lev[cols, l].astype(float)))
            self.grad.append(terms)

    def upto(self, N):
        """Exponents of all monomials with |g| <= N (the basis order used everywhere)."""
        return np.concatenate(self.levels[:N + 1])

    def evaluate(self, X, N):
        """Values of all monomials |g| <= N at the points X (m x n), ordered as upto(N)."""
        cols = [np.ones((len(X), 1))]
        for d in range(1, N + 1):
            cols.append(X[:, self.mul_axis[d]] * cols[-1][:, self.mul_pred[d]])
        return np.concatenate(cols, axis=1)

    def to_dense(self, per_level):
        """Level-wise moment vectors -> dense array mom[g] of shape (D+1,)*n."""
        mom = np.zeros((self.D + 1,) * self.n)
        for lev, vals in zip(self.levels, per_level):
            mom[tuple(lev.T)] = vals
        return mom


def ball_moments(n, D):
    """Exact moments int_{|x|<1} x^g dx of the unit n-ball (for the self-test)."""
    mom = np.zeros((D + 1,) * n)
    for g in np.ndindex(*mom.shape):
        if sum(g) <= D and all(e % 2 == 0 for e in g):
            mom[g] = math.exp(sum(math.lgamma((e + 1) / 2) for e in g)
                              - math.lgamma((sum(g) + n) / 2 + 1))
    return mom


# ─── polytope geometry & exact moments ────────────────────────────────────

CELL_NAMES = {
    ((3, 4),): "tetrahedron", ((3, 8),): "octahedron", ((4, 6),): "cube",
    ((3, 20),): "icosahedron", ((5, 12),): "dodecahedron",
    ((3, 8), (4, 6)): "cuboctahedron", ((3, 20), (5, 12)): "icosidodecahedron",
    ((3, 4), (6, 4)): "truncated tetrahedron", ((4, 6), (6, 8)): "truncated octahedron",
    ((3, 8), (8, 6)): "truncated cube", ((5, 12), (6, 20)): "truncated icosahedron",
    ((3, 20), (10, 12)): "truncated dodecahedron", ((3, 8), (4, 18)): "rhombicuboctahedron",
    ((3, 20), (4, 30), (5, 12)): "rhombicosidodecahedron",
    ((4, 12), (6, 8), (8, 6)): "truncated cuboctahedron",
    ((4, 30), (6, 20), (10, 12)): "truncated icosidodecahedron",
    ((3, 2), (4, 3)): "triangular prism", ((4, 5), (5, 2)): "pentagonal prism",
    ((4, 6), (6, 2)): "hexagonal prism", ((4, 8), (8, 2)): "octagonal prism",
    ((4, 10), (10, 2)): "decagonal prism", ((3, 8), (4, 2)): "square antiprism",
    ((3, 10), (5, 2)): "pentagonal antiprism", ((3, 4), (4, 1)): "square pyramid",
}


def _unit(v):
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def _plane_bases(points_list, dim):
    """Orthonormal bases (4 x dim) of the affine hulls of point sets (best-fit via SVD)."""
    out = [None] * len(points_list)
    by_size = {}
    for i, p in enumerate(points_list):
        by_size.setdefault(len(p), []).append(i)
    for idx in by_size.values():
        P = np.stack([points_list[i] for i in idx])
        P = P - P.mean(axis=1, keepdims=True)
        _, _, vt = np.linalg.svd(P, full_matrices=True)
        for k, i in enumerate(idx):
            out[i] = vt[k]              # rows: principal directions, last = normal(s)
    return out


class Polytope:
    """Boundary complex of a 4-polytope (topology_output JSON format)."""

    def __init__(self, path):
        with open(path, encoding="utf-8-sig") as f:
            d = json.load(f)
        self.path = path
        self.name = d["name"]
        self.description = d.get("description", "")
        X = np.asarray(d["vertices"], float)
        self.edges = np.asarray(d["edges"], dtype=np.intp)
        self.faces = [np.asarray(f, dtype=np.intp) for f in d["faces2d"]]
        self.cells = [np.asarray(c, dtype=np.intp) for c in d["cells"]]
        self.cell_faces = [np.asarray(c, dtype=np.intp) for c in d["cell_faces"]]
        file_normals = np.asarray(d["normals"], float)

        a, b = self.edges.T
        self.edge_length = float(np.linalg.norm(X[a] - X[b], axis=1).mean())
        Y = X / self.edge_length                          # physical units: edge = 1
        self.center = Y.mean(axis=0)
        self.radius = float(np.linalg.norm(Y - self.center, axis=1).max())
        P = self.xi = (Y - self.center) / self.radius     # Ritz coordinates, |xi| <= 1

        eid = {(min(u, v), max(u, v)): i for i, (u, v) in enumerate(self.edges)}
        self.face_edges = []
        for f in self.faces:
            fs = sorted(f.tolist())
            fe = [eid[(u, v)] for i, u in enumerate(fs) for v in fs[i + 1:] if (u, v) in eid]
            if len(fe) != len(f):
                raise ValueError(f"{self.name}: face {f.tolist()} has {len(fe)} edges")
            self.face_edges.append(np.asarray(fe, dtype=np.intp))

        # reference points (vertex centroids) and face / cell hyperplanes
        self.edge_mid = (P[a] + P[b]) / 2
        self.face_ctr = np.array([P[f].mean(axis=0) for f in self.faces])
        self.cell_ctr = np.array([P[c].mean(axis=0) for c in self.cells])
        face_basis = [vt[:2] for vt in _plane_bases([P[f] for f in self.faces], 2)]
        normals = np.array([vt[3] for vt in _plane_bases([P[c] for c in self.cells], 3)])
        normals *= np.sign(np.einsum("ij,ij->i", normals, file_normals))[:, None]
        self.cell_normal = normals
        self.h = np.einsum("ij,ij->i", normals, self.cell_ctr)   # signed facet distances

        # facet distances d(x0_G, H) for the Lasserre recursion
        self.edge_half = np.linalg.norm(P[b] - P[a], axis=1) / 2
        self.fe_dist = []
        for fi, fe in enumerate(self.face_edges):
            t = _unit(P[b[fe]] - P[a[fe]])
            w = self.edge_mid[fe] - self.face_ctr[fi]
            self.fe_dist.append(np.linalg.norm(w - np.einsum("ij,ij->i", w, t)[:, None] * t, axis=1))
        self.cf_dist = []
        for ci, cf in enumerate(self.cell_faces):
            w = self.face_ctr[cf] - self.cell_ctr[ci]
            Q = np.stack([face_basis[f] for f in cf])                     # (k, 2, 4)
            w_in = np.einsum("kpi,kp->ki", Q, np.einsum("kpi,ki->kp", Q, w))
            self.cf_dist.append(np.linalg.norm(w - w_in, axis=1))
        dmin = min(min(x.min() for x in self.fe_dist), min(x.min() for x in self.cf_dist))
        if dmin <= 0:
            raise ValueError(f"{self.name}: degenerate face/cell geometry")

    # -- symmetry -----------------------------------------------------------
    def symmetries(self, tol=1e-6):
        """Coordinate reflections x_l -> -x_l and the central inversion mapping P onto itself."""
        tree = cKDTree(self.xi)
        cellset = {frozenset(c.tolist()) for c in self.cells}

        def invariant(T):
            dist, perm = tree.query(self.xi * T)
            return dist.max() < tol and all(frozenset(perm[c].tolist()) in cellset
                                            for c in self.cells)

        refl = [l for l in range(4) if invariant(np.where(np.arange(4) == l, -1.0, 1.0))]
        return refl, invariant(-np.ones(4))

    # -- exact moments (Lasserre) --------------------------------------------
    def _cell_batches(self, max_rows):
        batch, rows = [], 0
        for c in range(len(self.cells)):
            r = 2 * len(self.cells[c]) + 2 * len(self.cell_faces[c]) - 1   # V+E+F+1 (Euler)
            if batch and rows + r > max_rows:
                yield np.array(batch)
                batch, rows = [], 0
            batch.append(c)
            rows += r
        if batch:
            yield np.array(batch)

    @staticmethod
    def _lift(facet_term, old, x0, grad_terms, k_plus_q):
        out = facet_term
        for l, cols, pred, coef in grad_terms:
            out[:, cols] += x0[:, l:l + 1] * (old[:, pred] * coef)
        out /= k_plus_q
        return out

    def _recursion(self, mono, cb):
        """Lasserre recursion for the cells cb: yields (d, Ic), Ic[i] = integrals of the
        degree-d monomials over cell cb[i] (3D measure)."""
        a, b = self.edges.T
        fb = np.unique(np.concatenate([self.cell_faces[c] for c in cb]))
        eb = np.unique(np.concatenate([self.face_edges[f] for f in fb]))
        vb = np.unique(self.edges[eb].ravel())
        ne, nf, nc = len(eb), len(fb), len(cb)
        rows = np.repeat(np.arange(ne), 2)
        cols = np.searchsorted(vb, np.stack([a[eb], b[eb]], axis=1).ravel())
        S_ev = sp.csr_matrix((np.repeat(self.edge_half[eb], 2), (rows, cols)), shape=(ne, len(vb)))
        rows = np.concatenate([np.full(len(self.face_edges[f]), i) for i, f in enumerate(fb)])
        cols = np.searchsorted(eb, np.concatenate([self.face_edges[f] for f in fb]))
        S_fe = sp.csr_matrix((np.concatenate([self.fe_dist[f] for f in fb]), (rows, cols)),
                             shape=(nf, ne))
        rows = np.concatenate([np.full(len(self.cell_faces[c]), i) for i, c in enumerate(cb)])
        cols = np.searchsorted(fb, np.concatenate([self.cell_faces[c] for c in cb]))
        S_cf = sp.csr_matrix((np.concatenate([self.cf_dist[c] for c in cb]), (rows, cols)),
                             shape=(nc, nf))
        Xv, Xe, Xf, Xc = self.xi[vb], self.edge_mid[eb], self.face_ctr[fb], self.cell_ctr[cb]

        Iv = np.ones((len(vb), 1))
        Ie = S_ev @ Iv
        If = (S_fe @ Ie) / 2
        Ic = (S_cf @ If) / 3
        yield 0, Ic
        for d in range(1, mono.D + 1):
            g = mono.grad[d]
            Iv = Xv[:, mono.mul_axis[d]] * Iv[:, mono.mul_pred[d]]
            Ie = self._lift(S_ev @ Iv, Ie, Xe, g, 1 + d)
            If = self._lift(S_fe @ Ie, If, Xf, g, 2 + d)
            Ic = self._lift(S_cf @ If, Ic, Xc, g, 3 + d)
            yield d, Ic

    def moments(self, mono, max_rows=4000):
        """int_P xi^g dxi for all |g| <= mono.D (dense array), P in Ritz coordinates."""
        acc = [np.zeros(len(lev)) for lev in mono.levels]
        for cb in self._cell_batches(max_rows):
            h = self.h[cb]
            for d, Ic in self._recursion(mono, cb):
                acc[d] += h @ Ic / (4 + d)       # x0 = origin: no gradient term
        return mono.to_dense(acc)

    def cell_moments(self, mono, c):
        """int_cell xi^g dS for one cell (dense array, Ritz coordinates, 3D measure)."""
        return mono.to_dense([Ic[0] for _, Ic in self._recursion(mono, np.array([c]))])

    # -- hit points -----------------------------------------------------------
    def vertex_normals(self):
        """Hit direction at a vertex: normalised sum of the incident cells' outward normals."""
        acc = np.zeros_like(self.xi)
        for c, verts in enumerate(self.cells):
            acc[verts] += self.cell_normal[c]
        return _unit(acc)

    def cell_label(self, c):
        census = {}
        for f in self.cell_faces[c]:
            census[len(self.faces[f])] = census.get(len(self.faces[f]), 0) + 1
        key = tuple(sorted(census.items()))
        return CELL_NAMES.get(key, "cell(" + ",".join(f"{k}^{m}" for k, m in key) + ")")


# ─── Rayleigh-Ritz ────────────────────────────────────────────────────────

def ritz_matrices(mom, E, lam, mu):
    """Scalar Gram matrix G and stiffness K (n*nb square) of the basis x^alpha e_i, alpha in E.

    a(u, v) = int lam div u div v + 2 mu eps(u):eps(v); for u = x^a e_i, v = x^b e_j:
    (lam a_i b_j + mu a_j b_i) m[a+b-e_i-e_j] + mu delta_ij sum_l a_l b_l m[a+b-2e_l].
    """
    nb, n = E.shape
    S = (E[:, None, :] + E[None, :, :]).astype(np.int32)

    def m(shift):
        idx = S - shift
        ok = (idx >= 0).all(axis=2)
        return np.where(ok, mom[tuple(np.maximum(idx, 0).transpose(2, 0, 1))], 0.0)

    eye = np.eye(n, dtype=np.intp)
    Ef = E.astype(float)
    G = m(np.zeros(n, dtype=np.intp))
    lap = sum(mu * np.outer(Ef[:, l], Ef[:, l]) * m(2 * eye[l]) for l in range(n))
    K = np.empty((n * nb, n * nb))
    for i in range(n):
        for j in range(n):
            blk = (lam * np.outer(Ef[:, i], Ef[:, j]) + mu * np.outer(Ef[:, j], Ef[:, i])) \
                  * m(eye[i] + eye[j])
            if i == j:
                blk += lap
            K[i * nb:(i + 1) * nb, j * nb:(j + 1) * nb] = blk
    return G, K


def ritz_solve(mom, mono, N, lam, mu, refl, central, count, vectors=True, rtol=1e-15):
    """Lowest `count` eigenpairs (rigid modes included) of K c = omega^2 M c.

    Returns omega^2 (ascending), coefficients C (n*nb x count, None unless vectors) of the
    M-orthonormal modes in the basis x^alpha e_i (component-major: row = i*nb + alpha),
    and the number of basis directions dropped as numerically dependent.
    """
    E = mono.upto(N)
    nb, n = E.shape
    G, K = ritz_matrices(mom, E, lam, mu)
    comp = np.repeat(np.arange(n), nb)
    mon = np.tile(np.arange(nb), n)
    alpha = E[mon]
    # parity of each basis field under the detected symmetries -> independent blocks
    keys = [(alpha[:, l] + (comp == l)) % 2 for l in refl]
    if central and len(refl) < n:
        keys.append((alpha.sum(axis=1) + 1) % 2)
    keys = np.stack(keys, axis=1) if keys else np.zeros((n * nb, 1), dtype=np.intp)
    _, block = np.unique(keys, axis=0, return_inverse=True)
    block = block.ravel()

    results, dropped = [], 0
    for bl in np.unique(block):
        idx = np.flatnonzero(block == bl)
        # canonical (Loewdin) orthogonalisation of the block-diagonal Gram matrix after
        # normalising every monomial (pentachoron, degree 12: raw cond ~1e17, ~3e13 after
        # scaling). Dropping directions would break the symmetry of the Ritz space and
        # split degenerate multiplets, so rtol only guards against true singularity.
        pieces, r = [], 0
        for i in range(n):
            loc = np.flatnonzero(comp[idx] == i)
            if len(loc) == 0:
                continue
            ms = mon[idx[loc]]
            g = G[np.ix_(ms, ms)]
            dsc = 1 / np.sqrt(np.diag(g))
            s, U = np.linalg.eigh(g * dsc[:, None] * dsc[None, :])
            keep = s > rtol * s[-1]
            pieces.append((loc, dsc[:, None] * U[:, keep] / np.sqrt(s[keep])))
            r += int(keep.sum())
            dropped += int((~keep).sum())
        X = np.zeros((len(idx), r))
        c0 = 0
        for loc, Xi in pieces:
            X[loc, c0:c0 + Xi.shape[1]] = Xi
            c0 += Xi.shape[1]
        A = X.T @ K[np.ix_(idx, idx)] @ X
        A = (A + A.T) / 2
        k = min(count, r)
        if vectors:
            w, Y = sla.eigh(A, subset_by_index=[0, k - 1])
            results.append((w, idx, X @ Y))
        else:
            results.append((sla.eigh(A, subset_by_index=[0, k - 1], eigvals_only=True), idx, None))

    w_all = np.concatenate([w for w, _, _ in results])
    src = np.concatenate([np.full(len(w), i) for i, (w, _, _) in enumerate(results)])
    col = np.concatenate([np.arange(len(w)) for w, _, _ in results])
    order = np.argsort(w_all, kind="stable")[:count]
    if not vectors:
        return w_all[order], None, dropped
    C = np.zeros((n * nb, len(order)))
    for j, o in enumerate(order):
        _, idx, Cb = results[src[o]]
        C[idx, j] = Cb[:, col[o]]
    return w_all[order], C, dropped


def group_multiplets(omega, rtol):
    """Split ascending frequencies into runs whose neighbours differ by < rtol (relative)."""
    groups, start = [], 0
    for k in range(1, len(omega) + 1):
        if k == len(omega) or omega[k] - omega[k - 1] > rtol * omega[k - 1]:
            groups.append(np.arange(start, k))
            start = k
    return groups


# ─── site gains ───────────────────────────────────────────────────────────

def site_gains(points, normals, mono, N, C, groups, chunk=4000):
    """Multiplet sums of squared mode displacements at points (Ritz-coordinate modes).

    gain_normal[s, g]     = sum_{k in g} (u_k(p_s) . n_s)^2
    gain_tangential[s, g] = sum_{k in g} |u_k(p_s) - (u_k . n_s) n_s|^2 / 3
    Both are invariant under the choice of basis inside a degenerate multiplet.
    """
    n = normals.shape[1]
    nb = C.shape[0] // n
    ind = np.zeros((C.shape[1], len(groups)))
    for g, members in enumerate(groups):
        ind[members, g] = 1.0
    gn = np.empty((len(points), len(groups)))
    gt = np.empty_like(gn)
    for s in range(0, len(points), chunk):
        Phi = mono.evaluate(points[s:s + chunk], N)
        U = [Phi @ C[i * nb:(i + 1) * nb] for i in range(n)]         # n x (m, modes)
        un = sum(U[i] * normals[s:s + chunk, i:i + 1] for i in range(n))
        uu = sum(Ui * Ui for Ui in U)
        gn[s:s + chunk] = (un * un) @ ind
        # difference of two sums: clamp the ~1e-16 rounding residue where it is exactly 0
        gt[s:s + chunk] = np.maximum(uu @ ind - gn[s:s + chunk], 0.0) / (n - 1)
    return gn, gt


def cell_mean_gains(poly, mono, N, C, groups, c):
    """Gains averaged over cell c (uniform hit point, direction = cell normal), exactly:
    int_cell (u.n)^2 dS = c_n^T G_cell c_n with the cell's monomial Gram matrix. Gain fields
    peak strongly at the cell's corners (up to ~150x the mean), so low-order quadrature fails;
    this agrees with a degree-25 quadrature to ~1e-11."""
    cm = poly.cell_moments(mono, c)
    E = mono.upto(N)
    nb = len(E)
    S = (E[:, None, :] + E[None, :, :]).astype(np.int32)
    Gc = cm[tuple(S.transpose(2, 0, 1))]
    area = cm[(0,) * 4]
    n = poly.cell_normal[c]
    Ci = [C[i * nb:(i + 1) * nb] for i in range(4)]
    Cn = sum(n[i] * Ci[i] for i in range(4))
    normal = np.einsum("ak,ak->k", Cn, Gc @ Cn) / area
    total = sum(np.einsum("ak,ak->k", X, Gc @ X) for X in Ci) / area
    ind = np.zeros((C.shape[1], len(groups)))
    for g, members in enumerate(groups):
        ind[members, g] = 1.0
    return normal @ ind, np.maximum(total - normal, 0.0) @ ind / 3, area


def classify(labels, feats, rtol):
    """Group sites with equal geometric label and (numerically) equal gain vectors."""
    scale = np.abs(feats).max() or 1.0
    cls = np.full(len(labels), -1)
    reps = []
    lab = np.asarray(labels, dtype=object)
    for s in range(len(labels)):
        if cls[s] >= 0:
            continue
        same = (lab == labels[s]) & (cls < 0)
        close = np.abs(feats - feats[s]).max(axis=1) <= rtol * scale
        cls[same & close] = len(reps)
        reps.append(s)
    return cls, reps


# ─── driver ───────────────────────────────────────────────────────────────

def analyse(poly, degree, lam, n_modes, max_error, mult_tol=1e-5, verbose=True):
    t0 = time.time()
    mono = Monomials(4, 2 * degree)
    mom = poly.moments(mono)
    refl, central = poly.symmetries()
    for l in refl:                                    # exact zeros instead of 1e-9 noise
        mom[(slice(None),) * l + (slice(1, None, 2),)] = 0.0
    if central:
        mom[np.indices(mom.shape).sum(axis=0) % 2 == 1] = 0.0
    t1 = time.time()

    count = n_modes + n_rigid(4)
    w, C, dropped = ritz_solve(mom, mono, degree, lam, 1.0, refl, central, count)
    w_lo, _, _ = ritz_solve(mom, mono, degree - 2, lam, 1.0, refl, central, count, vectors=False)
    if dropped:
        print(f"  WARNING {poly.name}: {dropped} numerically dependent basis directions dropped")
    t2 = time.time()

    nr = n_rigid(4)
    if not (np.abs(w[:nr]).max() < 1e-8 * w[nr]):
        raise RuntimeError(f"{poly.name}: rigid-body modes not separated: {w[:nr + 1]}")
    # the Ritz coordinates already have circumradius 1 (the output unit), c_s = 1
    omega = np.sqrt(w[nr:])
    omega_lo = np.sqrt(np.maximum(w_lo[nr:], 0))
    C = C[:, nr:]                                     # mass-normalised: int |u|^2 = 1
    rel_err = (omega_lo - omega) / omega              # >= 0 (nested Ritz spaces)

    # symmetry-degenerate modes agree to ~1e-9. Genuinely distinct multiplets closer than
    # mult_tol (e.g. 3.7e-6 in prico) are merged: inaudible, and their eigenvectors mix
    # numerically, so only the merged sums are exactly symmetric
    groups = group_multiplets(omega, mult_tol)
    if groups[-1][-1] == len(omega) - 1:              # possibly cut inside the multiplet
        groups = groups[:-1]
    keep = []
    for g in groups:
        if rel_err[g].max() > max_error:
            break
        keep.append(g)
    groups = keep
    C = C[:, : groups[-1][-1] + 1]

    volume = float(mom[(0,) * 4])
    mass = volume                                     # rho = 1

    # vertices: all equivalent in a uniform polychoron -> one gain table
    gn, gt = site_gains(poly.xi, poly.vertex_normals(), mono, degree, C, groups)
    feats = np.hstack([gn, gt])
    vspread = np.abs(feats - feats[0]).max() / np.abs(feats).max()
    vertex = (gn.mean(axis=0) * mass, gt.mean(axis=0) * mass)

    # cells: classes of symmetry-equivalent cells (same shape and same gains at the centre),
    # then the exact mean over one representative; a second member cross-checks the class
    labels = [poly.cell_label(c) for c in range(len(poly.cells))]
    gn, gt = site_gains(poly.cell_ctr, poly.cell_normal, mono, degree, C, groups)
    ccls, creps = classify(labels, np.hstack([gn, gt]), 1e-4)
    cells, check = [], vspread
    for k, r in enumerate(creps):
        members = np.flatnonzero(ccls == k)
        mn, mt, area = cell_mean_gains(poly, mono, degree, C, groups, r)
        if len(members) > 1:
            mn2, _, _ = cell_mean_gains(poly, mono, degree, C, groups, members[-1])
            check = max(check, np.abs(mn2 - mn).max() / mn.max())
        cells.append({"label": labels[r], "count": len(members), "area": area,
                      "normal": mn * mass, "tangential": mt * mass})
    # symmetry-equivalent vertices / cells must share their gains (distinct orbits differ ~1)
    if check > 1e-3:
        raise RuntimeError(f"{poly.name}: equivalent vertices or cells differ ({check:.1e})")
    if check > 1e-5:
        print(f"  WARNING {poly.name}: equivalent vertices or cells differ by {check:.1e}")
    same = {}
    for s in cells:
        same[s["label"]] = same.get(s["label"], 0) + 1
    seen = {}
    for s in cells:
        if same[s["label"]] > 1:
            seen[s["label"]] = seen.get(s["label"], 0) + 1
            s["label"] += f" #{seen[s['label']]}"
    wts = np.array([s["count"] * s["area"] for s in cells])
    mean = (sum(w * s["normal"] for w, s in zip(wts, cells)) / wts.sum(),
            sum(w * s["tangential"] for w, s in zip(wts, cells)) / wts.sum())
    t3 = time.time()

    freq = np.array([omega[g].mean() / (2 * math.pi) for g in groups])
    mult = [len(g) for g in groups]
    err = [float(rel_err[g].max()) for g in groups]
    if verbose:
        print(f"  {poly.name:9s} V={len(poly.xi):5d} C={len(poly.cells):4d}  "
              f"sym={''.join('xyzw'[l] for l in refl) or '-'}{'+inv' if central else ''}  "
              f"modes={sum(mult):4d} in {len(groups):3d} multiplets  "
              f"f1={freq[0]:.5f} (x{mult[0]})  f_max/f1={freq[-1] / freq[0]:.2f}  "
              f"cell classes={len(cells)} (check {check:.0e})  "
              f"[{t1 - t0:.1f}s moments, {t2 - t1:.1f}s eigen, {t3 - t2:.1f}s gains]", flush=True)
    return {
        "freq": freq, "mult": mult, "err": err, "volume": volume, "refl": refl, "mult_tol": mult_tol,
        "central": central, "mean": mean, "vertex": vertex, "cells": cells,
        "cell_class": ccls.tolist(),
    }


def sig(x, digits):
    return float(f"{x:.{digits}g}")


def gains(normal, tangential):
    return {"normal": [sig(x, 5) for x in normal], "tangential": [sig(x, 5) for x in tangential]}


def write_output(poly, res, out_dir, degree, lam, max_error):
    nb = math.comb(degree + 4, 4)
    doc = {
        "name": poly.name,
        "description": poly.description,
        "source": os.path.relpath(poly.path).replace("\\", "/"),
        "format": FORMAT,
        "model": {
            "body": "solid, homogeneous, isotropic, linear-elastic 4D body; free boundary",
            "equations": "rho u_tt = (lambda+mu) grad div u + mu Laplace u in R^4, traction-free",
            "method": f"Rayleigh-Ritz, complete polynomials of total degree {degree} per "
                      f"displacement component ({4 * nb} basis fields), exact polytope moments",
            "polynomial_degree": degree,
            "lambda_over_mu": lam,
            "units": "circumradius R = 1, density rho = 1, shear modulus mu = 1 (shear wave "
                     "speed c_s = 1, longitudinal c_p = sqrt(lambda/mu + 2))",
            "frequency": "f_Hz = frequency * c_s[m/s] / R[m]",
            "multiplet_merge_tol": res["mult_tol"],
            "gains": "mass * sum over the multiplet of (u_k(p) . d)^2 for mass-normalised modes "
                     "(int rho |u|^2 = 1): excitation by a unit impulse along d at p, picked up "
                     "along d, relative to the rigid-body value 1/mass. normal: d = outward "
                     "cell normal (vertex: normalised sum of the incident cell normals); "
                     "tangential: mean over the 3 tangent directions. cells[s]: p uniform over "
                     "a cell of class s; mean: p uniform over the whole boundary; vertex: p at "
                     "a vertex (all vertices are equivalent)",
            "max_rel_error": max_error,
            "rigid_body_modes_removed": n_rigid(4),
        },
        "geometry": {
            "volume": sig(res["volume"], 10),
            "edge_length": sig(1 / poly.radius, 10),
            "counts": [len(poly.xi), len(poly.edges), len(poly.faces), len(poly.cells)],
            "reflection_axes": [int(l) for l in res["refl"]],
            "centrally_symmetric": bool(res["central"]),
        },
        "modes": {
            "frequency": [sig(f, 8) for f in res["freq"]],
            "multiplicity": res["mult"],
            "rel_error": [sig(e, 2) for e in res["err"]],
        },
        "gains": {
            "mean": gains(*res["mean"]),
            "vertex": gains(*res["vertex"]),
            "cells": [dict(label=s["label"], count=s["count"],
                           **gains(s["normal"], s["tangential"])) for s in res["cells"]],
        },
        "cell_class": res["cell_class"],
    }
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, poly.name + ".json")
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(doc, f, separators=(",", ":"))
    os.replace(tmp, path)
    return path


# ─── self-test against analytic solutions ─────────────────────────────────

def _first_root(F, lo=0.5, hi=40.0, step=0.01):
    from scipy.optimize import brentq
    x = lo
    while x < hi:
        if F(x) * F(x + step) < 0:
            return brentq(F, x, x + step, xtol=1e-14)
        x += step
    raise ValueError("no root")


def selftest():
    from scipy.special import jv
    ok = True

    print("moments: tesseract (exact), pentachoron (cone quadrature)")
    for name in ("tes", "pen", "ico"):
        poly = Polytope(os.path.join("topology_output", name + ".json"))
        mono = Monomials(4, 12)
        mom = poly.moments(mono)
        if name == "tes":           # xi-cube [-1/2, 1/2]^4
            ref = np.zeros_like(mom)
            for g in np.ndindex(*mom.shape):
                if sum(g) <= 12 and all(e % 2 == 0 for e in g):
                    ref[g] = np.prod([0.5 ** e / (e + 1) for e in g])
        else:                       # flag simplices (0, cell, face, edge endpoints), Gauss
            t, wt = np.polynomial.legendre.leggauss(9)
            t, wt = (t + 1) / 2, wt / 2
            T = np.stack(np.meshgrid(t, t, t, t, indexing="ij"), -1).reshape(-1, 4)
            W = np.prod(np.stack(np.meshgrid(wt, wt, wt, wt, indexing="ij"), -1).reshape(-1, 4), 1)
            W *= T[:, 0] ** 3 * T[:, 1] ** 2 * T[:, 2]
            Z = np.cumprod(T, axis=1)                                # ordered simplex
            ref = np.zeros_like(mom)
            E = mono.upto(12)
            for c, cf in enumerate(poly.cell_faces):
                for f in cf:
                    for e in poly.face_edges[f]:
                        va, vb = poly.xi[poly.edges[e]]
                        V = np.stack([np.zeros(4), poly.cell_ctr[c], poly.face_ctr[f], va, vb])
                        D = V[1:] - V[:-1]
                        x = Z @ D
                        vals = mono.evaluate(x, 12).T @ W * abs(np.linalg.det(D))
                        ref[tuple(E.T)] += vals
        err = np.abs(mom - ref).max() / np.abs(ref).max()
        print(f"  {name}: volume {mom[0, 0, 0, 0] * poly.radius ** 4:.10f} (edge 1), "
              f"max rel. moment error {err:.1e}")
        ok &= err < 1e-12

    lam = mu = 1.0
    for n in (3, 4):
        N = 10
        mono = Monomials(n, 2 * N)
        mom = ball_moments(n, 2 * N)
        w, _, _ = ritz_solve(mom, mono, N, lam, mu, list(range(n)), True, 400, vectors=False)
        nr = n_rigid(n)
        omega = np.sqrt(np.maximum(w, 0))
        cs, cp = 1.0, math.sqrt((lam + 2 * mu))
        checks = []
        for l in (1, 2, 3):                                # torsional (toroidal) modes
            nu = (n - 2 + 2 * l) / 2
            z = _first_root(lambda z: (l - 1) * jv(nu, z) - z * jv(nu + 1, z))
            # multiplicity: 2l+1 in 3D; SO(4) irreps ((l+1)/2,(l-1)/2)+((l-1)/2,(l+1)/2) in 4D
            checks.append((f"torsional l={l}", z * cs, 2 * l + 1 if n == 3 else 2 * l * (l + 2)))
        nu = n / 2 - 1                                     # radial breathing mode
        z = _first_root(lambda z: (lam + 2 * mu) * z * jv(nu, z) - 2 * mu * (n - 1) * jv(nu + 1, z))
        checks.append(("breathing", z * cp, 1))
        print(f"{n}D unit ball, degree {N}: rigid modes {np.sum(w < 1e-9 * w[nr])} "
              f"(expected {nr}), lowest elastic omega {omega[nr]:.6f}")
        for label, ref, m_ref in checks:
            got = omega[np.argmin(np.abs(omega - ref))]
            mult = int(np.sum(np.abs(omega - got) < 1e-7 * got))
            e = abs(got - ref) / ref
            print(f"  {label:14s} analytic {ref:.8f}  ritz {got:.8f}  rel.err {e:.1e}  "
                  f"(x{mult}, expected x{m_ref})")
            ok &= e < 1e-4 and mult == m_ref
        ok &= int(np.sum(w < 1e-9 * w[nr])) == nr
    print("SELFTEST", "OK" if ok else "FAILED")
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("names", nargs="*", help="polychora to analyse (default: all)")
    ap.add_argument("--inputs", nargs="+", default=["topology_output"])
    ap.add_argument("--out", default="modal_output")
    ap.add_argument("--degree", type=int, default=12, help="polynomial degree N of the Ritz basis")
    ap.add_argument("--lame-ratio", type=float, default=1.0,
                    help="lambda/mu (1 = Cauchy solid, the isotropic central-force material)")
    ap.add_argument("--modes", type=int, default=600,
                    help="elastic modes (with multiplicity) to compute before truncation")
    ap.add_argument("--max-error", type=float, default=0.01,
                    help="export multiplets up to the first whose estimated relative frequency "
                         "error (degree N-2 vs N) exceeds this")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        sys.exit(0 if selftest() else 1)

    paths = []
    for d in args.inputs:
        paths += sorted(glob.glob(os.path.join(d, "*.json")))
    if args.names:
        by_name = {os.path.basename(p)[:-5]: p for p in paths}
        missing = [n for n in args.names if n not in by_name]
        if missing:
            sys.exit(f"unknown: {missing}")
        paths = [by_name[n] for n in args.names]

    print(f"modal analysis: degree {args.degree}, lambda/mu = {args.lame_ratio}, "
          f"{len(paths)} polychora -> {args.out}/")
    for p in paths:
        poly = Polytope(p)
        res = analyse(poly, args.degree, args.lame_ratio, args.modes, args.max_error)
        write_output(poly, res, args.out, args.degree, args.lame_ratio, args.max_error)


if __name__ == "__main__":
    main()

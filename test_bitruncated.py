import numpy as np
from scipy.spatial import ConvexHull
from collections import defaultdict, Counter
import itertools

sqrt2 = 2**0.5
b = 1.0 + sqrt2

def check_topology(pts_list, label=""):
    pts = list(set(tuple(round(x,8) for x in v) for v in pts_list))
    pts_arr = np.array(pts)
    try:
        hull = ConvexHull(pts_arr)
        eqs = hull.equations
        groups = defaultdict(set)
        for i, s in enumerate(hull.simplices):
            n = eqs[i,:4].copy()
            for j in range(4):
                if abs(n[j])>1e-8:
                    if n[j]<0: n=-n
                    break
            key=(tuple(np.round(n,8)), round(eqs[i,4],8))
            groups[key].update(s)
        C=len(groups); V=len(hull.vertices)
        sizes = Counter(len(v) for v in groups.values())
        print(f"{label}: V={V} C={C} sizes={sizes}")
        # compute E (edges)
        all_edges = set()
        for vl in groups.values():
            vl2 = sorted(vl); cp = pts_arr[vl2]
            try:
                h3 = ConvexHull(cp)
                fn = h3.equations[:,:3].copy()
                fg = defaultdict(set)
                for ii,s2 in enumerate(h3.simplices):
                    nc=fn[ii].copy(); nm=np.linalg.norm(nc)
                    if nm<1e-10: continue
                    nc/=nm
                    for jj in range(3):
                        if abs(nc[jj])>1e-8:
                            if nc[jj]<0: nc=-nc; break
                    fg[(tuple(np.round(nc,7)),round(h3.equations[ii,3],7))].update(vl2[v] for v in s2)
                for fv in fg.values():
                    fvl=list(fv)
                    if len(fvl)>=3:
                        fp = pts_arr[fvl]
                        try:
                            h2=ConvexHull(fp[:,:3])
                            for e in h2.simplices:
                                all_edges.add(tuple(sorted([fvl[e[0]],fvl[e[1]]])))
                        except: pass
            except: pass
        F_rough = sum(len(list(defaultdict(set,
            {(tuple(np.round(h3.equations[ii,:3]/np.linalg.norm(h3.equations[ii,:3]),7)),
              round(h3.equations[ii,3],7)): None}
        ).keys()) for vl in groups.values() for ii in range(1) if True) for _ in [None])
        print(f"       V={V} E={len(all_edges)} C={C}")
    except Exception as e:
        print(f"{label}: ERROR {e}")

# Bitruncated tesseract attempt: permutations of (+-1, +-1, +-1, +-(1+sqrt2))
# with even parity (product of signs > 0)
pts_bt = []
for b_pos in range(4):
    others = [i for i in range(4) if i != b_pos]
    for b_sign in [1, -1]:
        for s0, s1, s2 in itertools.product([-1,1], repeat=3):
            product = b_sign * s0 * s1 * s2
            if product > 0:  # even parity
                v = [0.0]*4
                v[b_pos] = b_sign * b
                v[others[0]] = s0
                v[others[1]] = s1
                v[others[2]] = s2
                pts_bt.append(tuple(round(x,8) for x in v))

pts_bt = list(set(pts_bt))
print(f"Vertex count: {len(pts_bt)}")
check_topology(pts_bt, "bitruncated_tesseract?")

# Also try the rectified tesseract for comparison
pts_rt = []
for zero_pos in range(4):
    others = [i for i in range(4) if i != zero_pos]
    for signs in itertools.product([-1,1], repeat=3):
        v = [0.0]*4
        v[zero_pos] = 0
        for k,s in zip(others, signs):
            v[k] = s * sqrt2
        pts_rt.append(tuple(round(x,8) for x in v))
pts_rt = list(set(pts_rt))
print(f"\nRectified tesseract vertex count: {len(pts_rt)}")
check_topology(pts_rt, "rectified_tesseract")

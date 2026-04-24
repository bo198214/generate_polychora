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
            key=(tuple(np.round(n,7)), round(eqs[i,4],7))
            groups[key].update(s)
        C=len(groups); V=len(hull.vertices)
        sizes = Counter(len(v) for v in groups.values())
        print(f"{label}: V={V} C={C} sizes={sizes}")
        return V, C, pts
    except Exception as e:
        print(f"{label}: ERROR {e} V={len(pts)}")
        return len(pts), 0, pts

# Attempt 1: permutations of (0, +/-1, +/-1, +/-(1+sqrt2)) - 4*3*4*2=96 no
# Try with zero in specific places...

# The bitruncated tesseract in standard coordinates uses:
# permutations of (0, +/-1, +/-1, +/-(sqrt2+1)) but that's 96
#
# From CDy analysis: bitruncated = nodes 1,2 active in B4.
# Let me compute it using our PolychoraGenerator with nodes 1,2 active (an=6):
# That gives V=96 (wrong).
#
# Alternative: maybe the bitruncated tesseract vertices are at
# all permutations of (+/-a, +/-a, +/-b, +/-b) with a=1, b=1+sqrt2 AND
# with the constraint a*a + b*b = some normalization

# Actually: maybe the bitruncated tesseract has a very different set.
# From Coxeter: bitruncated tesseract = 2t{4,3,3}
# Vertices: permutations of (0, +/-1, +/-(1+sqrt2), +/-(1+sqrt2))
# Count: 4 * 2 * 4 * 2... let me compute

pts1 = []
for zero_pos in range(4):
    for s_a in [-1, 1]:
        for s_b1 in [-1, 1]:
            for s_b2 in [-1, 1]:
                remaining = [i for i in range(4) if i != zero_pos]
                v = [0.0]*4
                v[zero_pos] = 0
                # Pick one position for +/-1, two for +/-(1+sqrt2)
                for a_pos_idx in range(3):
                    a_pos = remaining[a_pos_idx]
                    b_positions = [remaining[j] for j in range(3) if j != a_pos_idx]
                    v2 = v.copy()
                    v2[a_pos] = s_a * 1.0
                    v2[b_positions[0]] = s_b1 * b
                    v2[b_positions[1]] = s_b2 * b
                    pts1.append(tuple(round(x,8) for x in v2))

check_topology(pts1, "zero+1+b+b")

# Attempt: (0, +/-1, +/-(1+sqrt2), +/-(1+sqrt2)) only two b's
# That's what I just computed - should give V = 4*3*2*4 = 96... let me check unique
pts2 = []
for perm in itertools.permutations([0, 1, b, b]):
    for signs in itertools.product([-1,1,1,1,1], repeat=1):
        pass  # too complex

# Try another approach: compute using Wythoff formula
# Bitruncated tesseract generator: nodes 0,3 INACTIVE, nodes 1,2 ACTIVE in B4
# This is the "2t" polytope
# From our formula: an=6 gives generator P=(0, sqrt2, sqrt2*(1+sqrt2)/?, ?)

# What if bitruncated tesseract vertices are permutations of (+/-1, +/-1, 0, 0) scaled differently?
# Or: (+/-c, +/-c, +/-c, +/-(1+sqrt2)) for some c?
# Count: 4*2^4 = 64 - that's not 32

# Let me try the ACTUAL coordinates from the Wythoff construction in a different way
# Using B4 matrix but with different Coxeter diagram orientation

# MAYBE: the correct B4 matrix for bitruncated should be
# {1, 3, 3, 4} type? Let me just try generating with H4 family looking for V=32
from scipy.spatial import ConvexHull

# Import what PolychoraGenerator does
def poly_gen(mat, an):
    n = len(mat)
    normals = [[0.0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i):
            dot = -np.cos(np.pi / mat[i][j])
            s = sum(normals[i][k]*normals[j][k] for k in range(j))
            normals[i][j] = (dot-s)/(normals[j][j] if normals[j][j]!=0 else 1)
        ssq = sum(normals[i][k]**2 for k in range(i))
        normals[i][i] = max(0,1.0-ssq)**0.5

    d = [(1.0 if (an>>i)&1 else 0.0) for i in range(n)]
    p = [0.0]*n
    for i in range(n):
        s = sum(normals[i][j]*p[j] for j in range(i))
        p[i] = (d[i]-s)/(normals[i][i] if normals[i][i]!=0 else 1)

    verts=[np.array(p)]; seen={tuple(np.round(p,8)+0.0)}
    i=0
    while i<len(verts):
        v=verts[i]; i+=1
        for mi in normals:
            dot=sum(v[k]*mi[k] for k in range(4))
            if abs(dot)<1e-8: continue
            nxt=v-2*dot*np.array(mi)
            key=tuple((np.round(nxt,8)+0.0).tolist())
            if key not in seen: seen.add(key); verts.append(nxt)
    return verts

# Try F4 family looking for V=32
f4mat = [[1,3,2,2],[3,1,4,2],[2,4,1,3],[2,2,3,1]]
for an in range(1,16):
    v = poly_gen(f4mat, an)
    if len(v)==32:
        print(f"F4 an={an}: V={len(v)}, first={list(np.round(v[0],4))}")
        check_topology(v, f"F4/an={an}")

import numpy as np
from scipy.spatial import ConvexHull
from collections import defaultdict

def get_mirror_normals(matrix):
    n = len(matrix)
    normals = [[0.0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i):
            dot = -np.cos(np.pi / matrix[i][j])
            s = sum(normals[i][k]*normals[j][k] for k in range(j))
            normals[i][j] = (dot - s) / (normals[j][j] if normals[j][j] != 0 else 1)
        ssq = sum(normals[i][k]**2 for k in range(i))
        normals[i][i] = max(0.0, 1.0 - ssq)**0.5
    return normals

def solve_gen(mirrors, active):
    n = len(mirrors)
    d = [(1.0 if (active>>i)&1 else 0.0) for i in range(n)]
    p = [0.0]*n
    for i in range(n):
        s = sum(mirrors[i][j]*p[j] for j in range(i))
        p[i] = (d[i]-s)/(mirrors[i][i] if mirrors[i][i]!=0 else 1)
    return p

def gen_orbit(start, mirrors):
    verts=[np.array(start)]; seen={tuple(np.round(start,8))}
    i=0
    while i<len(verts):
        v=verts[i]; i+=1
        for mi in mirrors:
            dot=sum(v[k]*mi[k] for k in range(4))
            if abs(dot)<1e-8: continue
            nxt=v - 2*dot*np.array(mi)
            key=tuple(np.round(nxt,8))
            if key not in seen: seen.add(key); verts.append(nxt)
    return np.array(verts)

def hull_topology(pts):
    try:
        h=ConvexHull(pts)
        eqs=h.equations
        groups=defaultdict(set)
        for i,s in enumerate(h.simplices):
            n=eqs[i,:4].copy()
            for j in range(4):
                if abs(n[j])>1e-8:
                    if n[j]<0: n=-n
                    break
            key=(tuple(np.round(n,7)),round(eqs[i,4],7))
            groups[key].update(s)
        C=len(groups)
        V=len(h.vertices)
        all_faces=set()
        all_edges=set()
        for vl in groups.values():
            vl=sorted(vl); cp=pts[vl]
            try:
                h3=ConvexHull(cp)
                fn=h3.equations[:,:3].copy()
                fg=defaultdict(set)
                for i2,s2 in enumerate(h3.simplices):
                    nc=fn[i2].copy(); nm=np.linalg.norm(nc)
                    if nm<1e-10: continue
                    nc/=nm
                    for j in range(3):
                        if abs(nc[j])>1e-8:
                            if nc[j]<0: nc=-nc
                            break
                    fg[(tuple(np.round(nc,6)),round(h3.equations[i2,3],6))].update(vl[v] for v in s2)
                for fv in fg.values():
                    fk=tuple(sorted(fv)); all_faces.add(fk)
                    fvl=list(fv)
                    if len(fvl)>=3:
                        fp=pts[fvl]
                        try:
                            h2=ConvexHull(fp[:,:3])
                            for e in h2.simplices:
                                all_edges.add(tuple(sorted([fvl[e[0]],fvl[e[1]]])))
                        except: pass
            except: pass
        return V, len(all_edges), len(all_faces), C
    except: return 0,0,0,0

groups = {
    'A4': [[1,3,2,2],[3,1,3,2],[2,3,1,3],[2,2,3,1]],
    'B4': [[1,4,2,2],[4,1,3,2],[2,3,1,3],[2,2,3,1]],
    'F4': [[1,3,2,2],[3,1,4,2],[2,4,1,3],[2,2,3,1]],
    'H4': [[1,5,2,2],[5,1,3,2],[2,3,1,3],[2,2,3,1]],
}

expected = {
    'tes':(16,32,24,8),'hex':(8,24,32,16),'tah':(64,128,80,16),
    'rico':(32,96,80,16),'tico':(48,120,88,16),
    'sadi':(20,60,60,20),'ico':(24,96,96,24),'hi':(600,1200,720,120),
}

for gname, mat in groups.items():
    mir = get_mirror_normals(mat)
    print(f'\n{gname} family:')
    for an in range(1,16):
        p = solve_gen(mir, an)
        pts = gen_orbit(p, mir)
        vefc = hull_topology(pts)
        match = [k for k,v in expected.items() if v==vefc]
        tag = '  <- ' + ','.join(match) if match else ''
        print(f'  an={an:2d} V={vefc[0]:5d} E={vefc[1]:5d} F={vefc[2]:5d} C={vefc[3]:4d}{tag}')
    if gname == 'B4':
        break  # limit output for now

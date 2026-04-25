#!/usr/bin/env python3
"""
Integrity check for topology_output/ JSON files.
Usage: python check_topology.py [topology_dir]

Checks:
  1. Required fields present
  2. Euler characteristic V - E + F - C = 0
  3. normals count == cells count
  4. cell_faces length == cells count; all indices in range
  5. Closed-manifold condition via cell_faces: every ridge (2-face) is shared by
     exactly 2 cells (necessary & sufficient for a closed 3-manifold without
     boundary — stronger than Euler alone).
  6. Vertex consistency: face vertices are a subset of the owning cell's vertices.
"""
import json, os, sys
from collections import defaultdict

tdir = sys.argv[1] if len(sys.argv) > 1 else "topology_output"
required = ["name", "description", "vertices", "edges", "faces2d", "cells", "normals", "cell_faces"]

ok = bad = 0
for fname in sorted(os.listdir(tdir)):
    if not fname.endswith(".json"):
        continue
    path = os.path.join(tdir, fname)
    name = fname[:-5]
    try:
        with open(path, encoding="utf-8-sig") as f:
            d = json.load(f)

        missing = [k for k in required if k not in d]
        v  = len(d.get("vertices",  []))
        e  = len(d.get("edges",     []))
        fa = len(d.get("faces2d",   []))
        c  = len(d.get("cells",     []))
        n  = len(d.get("normals",   []))
        euler = v - e + fa - c

        issues = []
        if missing:
            issues.append(f"missing fields: {missing}")
        if euler != 0:
            issues.append(f"Euler={euler} != 0")
        if n != c:
            issues.append(f"normals={n} != cells={c}")

        if "cell_faces" in d:
            cf = d["cell_faces"]
            if len(cf) != c:
                issues.append(f"cell_faces length={len(cf)} != cells={c}")
            else:
                # All face indices in range
                bad_idx = [(ci, fi) for ci, faces in enumerate(cf)
                           for fi in faces if not (0 <= fi < fa)]
                if bad_idx:
                    issues.append(f"cell_faces out-of-range indices: {bad_idx[:5]}")

                # Each face must appear in exactly 2 cells (closed-manifold via cell_faces).
                # Use a full-size array so faces with 0 appearances are also caught.
                face_cell_count = [0] * fa
                for faces in cf:
                    for fi in faces:
                        if 0 <= fi < fa:
                            face_cell_count[fi] += 1
                bad_faces = [(fi, face_cell_count[fi]) for fi in range(fa)
                             if face_cell_count[fi] != 2]
                if bad_faces:
                    counts = defaultdict(int)
                    for _, cnt in bad_faces:
                        counts[cnt] += 1
                    details = ", ".join(f"{n} faces in {k} cells"
                                       for k, n in sorted(counts.items()))
                    issues.append(f"non-manifold (cell_faces): {details}")

                # Vertex consistency: face vertices must be subset of cell vertices
                cells_vtx = [set(cell) for cell in d["cells"]]
                faces_vtx = [set(face) for face in d["faces2d"]]
                incons = [(ci, fi) for ci, faces in enumerate(cf)
                          for fi in faces if not faces_vtx[fi].issubset(cells_vtx[ci])]
                if incons:
                    issues.append(f"cell_faces vertex inconsistency: {incons[:5]}")

        if issues:
            print(f"  FAIL {name:10}  V={v:5} E={e:5} F={fa:5} C={c:4}  {', '.join(issues)}")
            bad += 1
        else:
            desc = d.get("description", "")
            print(f"  ok   {name:10}  V={v:5} E={e:5} F={fa:5} C={c:4}  {desc[:50]}")
            ok += 1

    except Exception as ex:
        print(f"  ERROR {name}: {ex}")
        bad += 1

print(f"\n{ok} ok, {bad} failed")
sys.exit(1 if bad else 0)

#!/usr/bin/env python3
"""
Integrity check for topology_output/ JSON files.
Usage: python check_topology.py [topology_dir]

Checks:
  1. Required fields present
  2. Euler characteristic V - E + F - C = 0
  3. normals count == cells count
  4. Closed-manifold condition: every ridge (2-face) is shared by
     exactly 2 cells (necessary & sufficient for a closed 3-manifold
     without boundary — stronger than Euler alone).
"""
import json, os, sys
from collections import defaultdict

tdir = sys.argv[1] if len(sys.argv) > 1 else "topology_output"
required = ["name", "description", "vertices", "edges", "faces2d", "cells", "normals"]

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

        # Closed-manifold check: every ridge must belong to exactly 2 cells.
        # Build vertex→cells index, then intersect for each ridge.  O(C·|cell| + F·|face|).
        if not missing and c > 0 and fa > 0:
            vtx_to_cells = defaultdict(set)
            for ci, cell in enumerate(d["cells"]):
                for vtx in cell:
                    vtx_to_cells[vtx].add(ci)

            ridge_counts = defaultdict(int)  # count of cells per ridge
            for face in d["faces2d"]:
                if not face:
                    continue
                # Cells that contain ALL vertices of this ridge
                containing = vtx_to_cells[face[0]].copy()
                for vtx in face[1:]:
                    containing &= vtx_to_cells[vtx]
                ridge_counts[len(containing)] += 1

            non_two = {k: v for k, v in ridge_counts.items() if k != 2}
            if non_two:
                details = ", ".join(f"{cnt} ridges in {k} cells" for k, cnt in sorted(non_two.items()))
                issues.append(f"non-manifold: {details}")

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

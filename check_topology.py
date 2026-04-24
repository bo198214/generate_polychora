#!/usr/bin/env python3
"""
Integrity check for topology_output/ JSON files.
Usage: python check_topology.py [topology_dir]
"""
import json, os, sys

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

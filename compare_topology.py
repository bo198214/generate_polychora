"""
Compare topology of two polychoron JSON files.

Strategy: canonicalize both JSONs by re-indexing vertices in coordinate order
(rounded), then compare the resulting sorted edge/face/cell sets.
Falls back to a structural summary if coordinates differ too much.

Usage: python compare_topology.py a.json b.json
"""

import json
import sys
import math
from collections import Counter

COORD_PRECISION = 5  # decimal places for rounding


def load(path):
    with open(path) as f:
        return json.load(f)


def round_coord(c):
    return tuple(round(x, COORD_PRECISION) for x in c)


def canonical_index_map(vertices):
    """Sort vertices lexicographically by rounded coords → new index order."""
    indexed = sorted(enumerate(vertices), key=lambda iv: round_coord(iv[1]))
    # old_idx → new_idx
    return {old: new for new, (old, _) in enumerate(indexed)}


def reindex(index_list, imap):
    return tuple(sorted(imap[v] for v in index_list))


def canonical_sets(data):
    imap = canonical_index_map(data["vertices"])
    edges  = frozenset(reindex(e, imap) for e in data.get("edges",  []))
    faces  = frozenset(reindex(f, imap) for f in data.get("faces2d", []))
    cells  = frozenset(reindex(c, imap) for c in data.get("cells",  []))
    return edges, faces, cells


def structural_summary(data):
    """Counts + degree sequences — coordinate-independent."""
    n_v = len(data.get("vertices", []))
    n_e = len(data.get("edges",    []))
    n_f = len(data.get("faces2d",  []))
    n_c = len(data.get("cells",    []))

    vertex_edge_degree  = Counter(v for e in data.get("edges",   []) for v in e)
    vertex_face_degree  = Counter(v for f in data.get("faces2d", []) for v in f)
    vertex_cell_degree  = Counter(v for c in data.get("cells",   []) for v in c)
    face_sizes          = Counter(len(f) for f in data.get("faces2d", []))
    cell_sizes          = Counter(len(c) for c in data.get("cells",   []))

    return {
        "counts": (n_v, n_e, n_f, n_c),
        "edge_degree_seq":  tuple(sorted(vertex_edge_degree.values())),
        "face_degree_seq":  tuple(sorted(vertex_face_degree.values())),
        "cell_degree_seq":  tuple(sorted(vertex_cell_degree.values())),
        "face_sizes":       dict(sorted(face_sizes.items())),
        "cell_sizes":       dict(sorted(cell_sizes.items())),
    }


def compare(path_a, path_b):
    a = load(path_a)
    b = load(path_b)

    print(f"A: {path_a}  ({a.get('name','?')}) — {a.get('description','')}")
    print(f"B: {path_b}  ({b.get('name','?')}) — {b.get('description','')}")
    print()

    # --- Structural summary (always) ---
    sa = structural_summary(a)
    sb = structural_summary(b)

    labels = ["vertices", "edges", "faces2d", "cells"]
    print("Counts:")
    for label, ca, cb in zip(labels, sa["counts"], sb["counts"]):
        mark = "OK" if ca == cb else "DIFF"
        print(f"  {label:12s}  A={ca:5d}  B={cb:5d}  [{mark}]")
    print()

    for key in ("edge_degree_seq", "face_degree_seq", "cell_degree_seq",
                "face_sizes", "cell_sizes"):
        match = sa[key] == sb[key]
        print(f"  {key:24s}  [{'OK' if match else 'DIFF'}]")
        if not match:
            print(f"    A: {sa[key]}")
            print(f"    B: {sb[key]}")
    print()

    structural_ok = (sa == sb)
    print(f"Structural summary: {'IDENTICAL' if structural_ok else 'DIFFERENT'}")

    if not structural_ok:
        print("Topologies differ structurally — skipping canonical comparison.")
        return

    # --- Canonical comparison (requires matching coordinates) ---
    try:
        ea, fa, ca_ = canonical_sets(a)
        eb, fb, cb_ = canonical_sets(b)
    except Exception as ex:
        print(f"Canonical comparison failed: {ex}")
        return

    edges_ok = (ea == eb)
    faces_ok = (fa == fb)
    cells_ok = (ca_ == cb_)
    print(f"\nCanonical (coordinate-based) comparison:")
    print(f"  edges:   [{'OK' if edges_ok else 'DIFF'}]")
    print(f"  faces2d: [{'OK' if faces_ok else 'DIFF'}]")
    print(f"  cells:   [{'OK' if cells_ok else 'DIFF'}]")

    if not edges_ok:
        only_a = ea - eb
        only_b = eb - ea
        print(f"    only in A: {len(only_a)} edges, e.g. {list(only_a)[:3]}")
        print(f"    only in B: {len(only_b)} edges, e.g. {list(only_b)[:3]}")
    if not faces_ok:
        only_a = fa - fb
        only_b = fb - fa
        print(f"    only in A: {len(only_a)} faces, e.g. {list(only_a)[:3]}")
        print(f"    only in B: {len(only_b)} faces, e.g. {list(only_b)[:3]}")
    if not cells_ok:
        only_a = ca_ - cb_
        only_b = cb_ - ca_
        print(f"    only in A: {len(only_a)} cells, e.g. {list(only_a)[:3]}")
        print(f"    only in B: {len(only_b)} cells, e.g. {list(only_b)[:3]}")

    if edges_ok and faces_ok and cells_ok:
        print("\nResult: TOPOLOGIES ARE IDENTICAL")
    else:
        print("\nResult: TOPOLOGIES DIFFER")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python compare_topology.py a.json b.json")
        sys.exit(1)
    compare(sys.argv[1], sys.argv[2])

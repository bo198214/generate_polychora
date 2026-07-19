# generate_polychora

Generates the 47 non-prismatic convex uniform 4-polytopes as JSON (vertices + full
boundary topology) for the Tesserian Unity project. Two stages:

1. **Vertex generation** (`dotnet run -- vertices`): Wythoff construction — mirrors from
   the Coxeter matrix (A4/B4/F4/H4), generator point solved per active-node bitmask,
   vertex orbit by reflection closure (`PolychoraGenerator.cs`). The two non-Wythoffian
   polychora (sadi = snub 24-cell, gap = grand antiprism) come from `SnubGenerator.cs`.
   Output: `vertex_output/<name>.json`.
2. **Topology** (`dotnet run -- topology` / `make topology`): 4D gift-wrapping convex hull
   (`TrueConvexHull4D.cs`) computes cells, faces, edges, outward normals, cell→face
   incidence. Output: `topology_output/<name>.json`.

The finished `topology_output/*.json` files are copied into the Unity repo at
`tesserian/Assets/_Tesserian/RotatingPolychoron/Resources/polychora/` (the Polychoron Watch
loads all of them; content is kept byte-identical apart from the `name`/`description`
header fields).

Tests: `dotnet test Tests/Tests.csproj` verifies V/E/F/C of every generated file against
the literature values for all 47 polytopes (Klitzing / Wikipedia "Uniform 4-polytope").

## Naming: Bowers acronyms

File names and JSON `name` fields are the standard **Bowers acronyms** exactly as listed on
[Wikipedia's Uniform 4-polytope article](https://en.wikipedia.org/wiki/Uniform_4-polytope)
and [Klitzing's site](https://bendwavy.org/klitzing/dimensions/polychora.htm).

History (2026-07-19): the original name table used a mix of invented acronyms (`tappy`,
`dappat`, `scic`, `thic`, `xic`, `drico`, `rhi`, `rex`) and real acronyms that belong to
*different* polytopes (`hap` = hexagonal antiprism, `snic` = snub cube, `gic` = great
tetracontoctachoron, `frico` = facetorectified icositetrachoron, `trico` = triangle–24-cell
duoprism; `rico`/`tico`/`cont`/`tah`/`spic` were valid acronyms sitting on the wrong
neighbours within the 47). Every file's true identity was verified from its element counts
(V/E/F/C uniquely identify each of the 47) plus face-polygon census for the two count-ties
(prit vs proh via octagons, prix vs prahi via decagons), then renamed in both repos.
The descriptions in the old table were almost all correct — only the acronyms were off —
which is what made the mislabeling stable for so long.

## The prahi bug (pivot sweep orientation)

`prahi` (runcitruncated 600-cell) was missing for a long time, blamed on "numerical
instability". The real cause was **not precision** but a logic bug in the gift-wrapping
pivot:

- The pivot rotates a supporting hyperplane around a ridge (2-face) from the known cell to
  the adjacent one, picking the candidate with the smallest positive rotation angle.
- The rotation-plane basis vector was `p2 = cross4(e1, e2, prevNormal)` where `e1`, `e2`
  come from the first three ridge vertices **in arbitrary order** — so the sweep direction
  was effectively a coin flip per ridge. With an inverted sweep the true next cell sits at
  2π−θ, hidden behind thousands of interior planes; such ridges always failed.
- The algorithm still worked for almost everything because every cell has many inbound
  ridges and only needs one correctly-oriented ridge to be discovered. Exactly one
  hexagonal prism of prahi had **all 8** inbound ridges inverted → 2639 of 2640 cells,
  with all 13440 faces present and 8 of them dangling (incident to only one cell).
- Diagnosis was conclusive because the missing cell's hyperplane had margins of 0.22 vs
  noise of 1e-9 — a tolerance problem was impossible; arbitrary-precision arithmetic would
  not have helped.

Fix: each queued ridge carries a **reference vertex** of the discovering cell (off the
ridge), and `Pivot` orients `p2` so that reference point lies on the correct side
(`p2·(ref−v0) > 0`). The sign test is exact in the relevant regime (the reference point
lies *on* the previous cell's hyperplane, so the dot product is a pure pencil-plane
component of magnitude ~1). Additionally the pivot now walks the angle-sorted candidate
list and takes the first true supporting hyperplane instead of giving up when the very
smallest candidate is a non-extreme noise plane, and candidates coplanar with the previous
cell are excluded by direction (`cosθ > 1−1e-7`) rather than by signed angle, which
floating noise can push to θ ≈ +1e-9.

**Fail fast:** `TrueConvexHull4D.Compute` now throws if any face is not shared by exactly
two cells — the invariant that silently broke for prahi. For diagnosis of a broken hull,
set the environment variable `HULL_KEEP_BROKEN=1` to write the file anyway.

After the fix, all 47 topologies validate: literature element counts, Euler characteristic
0, every face shared by exactly two cells. Regenerating previously-good polytopes with the
fixed pivot yields set-identical topology (cell/face discovery order may differ).

## Nonconvex regular-faced polychora (excavation)

`dotnet run -- excavate` (`Excavation.cs`) applies **cell excavation** — the elementary
nonconvex CSG step: one boundary cell is replaced by the lateral cells of a unit-edged
pyramid whose apex points into the solid. This is the 4D version of Bonnie Stewart's
excavation move (the operation behind the 3D Stewart toroids) and the mirror image of the
"augmentation" used by the CRF community (hi.gher.space; see qfbox.info/4d/crf). The
convex CRF world is systematically explored (all regular/uniform polychora known, 314
million non-adjacent 600-cell diminishings enumerated); embedded *nonconvex* regular-faced
polychora are essentially uncharted — the community's nonconvex work concentrates on
self-intersecting star polychora instead.

Feasibility rules (checked at runtime): the pyramid needs base circumradius < edge
(icosahedron 0.951 ✓, octahedron 0.707 ✓, tetrahedron 0.612 ✓, dodecahedron 1.40 ✗ —
dodecahedral cells cannot be excavated unit-edged), and the apex must stay strictly inside
every original cell hyperplane (rules out the 16-cell: dent depth 1.118e > thickness 1.0e).
For a convex source the dent pyramid is automatically contained in the solid, so the
result is embedded; closedness and Euler characteristic are re-validated after surgery.

Curated outputs (`crf_output/`, copied to the Unity repo's
`Assets/_Tesserian/RotatingPolychoron/Resources/complexes/`, shown by the Watch's
dev Complex mode which renders nonconvex figures with real occlusion):

| file | source | surgery | V/E/F/C |
|---|---|---|---|
| excavated-ex | 600-cell | one tet cell → pentachoron dimple | 121/724/1206/603 |
| bi-excavated-ex | 600-cell | two antipodal pentachoron dimples | 122/728/1212/606 |
| excavated-ico | 24-cell | one octahedral-pyramid dent, apex exactly at the center | 25/102/108/31 |
| excavated-sadi | snub 24-cell | one shallow icosahedral-pyramid dimple (h ≈ 0.31e) | 97/444/510/163 |

Next steps if desired: multi-excavations with pairwise-disjoint pyramids, tunnel/toroid
constructions (boundary topology beyond the 3-sphere), and general boolean CSG.

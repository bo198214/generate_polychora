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
3. **Modal analysis** (`make modal`): vibration modes of every polychoron as a solid
   elastic 4D body (`modal_analysis.py`), plus example sounds (`make sounds`,
   `synth_modal.py`). Output: `modal_output/<name>.json`, `modal_output/sounds/*.wav`.
   See [Modal analysis](#modal-analysis-4d-vibration-modes).
4. **Hollow bodies** (`make hollow`): the same for thin-walled hollow polychora
   (`modal_hollow.py`, uses the cell classes from step 3). Output:
   `modal_hollow_output/<name>.json`. See [Hollow bodies](#hollow-bodies-thin-walled).

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

## Modal analysis (4D vibration modes)

`python modal_analysis.py` (or `make modal`, which runs `--selftest` first) computes the free
vibrations of every polychoron in `topology_output/`, treated as a **solid, homogeneous,
isotropic, linear-elastic 4D body** that is unsupported in R^4. Displacements u(x,t) ∈ R^4
obey the 4D Navier–Cauchy equations

    rho u_tt = (lambda + mu) grad(div u) + mu Laplace(u)    inside the polytope,
    sigma(u) n = 0                                          on its 3D boundary (the cells),

with sigma = lambda tr(eps) I + 2 mu eps. Default material: lambda = mu (the Cauchy solid,
`--lame-ratio`). Units: **circumradius R = 1** (the Polychoron Watch shows every polychoron
at radius 1, while `topology_output` has edge length 2), density 1, shear modulus 1, so the
shear wave speed is c_s = 1 and all frequencies are dimensionless:
**f_Hz = frequency · c_s / R**. The 10 rigid-body modes of a free 4D body (4 translations +
6 rotation planes) are removed.

**Method.** Rayleigh–Ritz with a complete polynomial basis — every displacement component
is a polynomial of total degree ≤ 12 in (x, y, z, w), 4 × 1820 basis fields; this is the 4D
version of Visscher's "xyz algorithm" from resonant ultrasound spectroscopy. No mesh: mass
and stiffness entries are integrals of monomials over the polytope, computed exactly with
Lasserre's divergence-theorem recursion vertex → edge → face → cell → polytope
((k+q)∫_G f = Σ_H dist(x0,H) ∫_H f + ∫_G x0·∇f for f homogeneous of degree q). Coordinate
reflections and central inversion of the polytope split the eigenproblem into parity
blocks. Two numerical details matter: every monomial is normalised before the Löwdin
orthogonalisation (degree 12 on the thin 5-cell has Gram condition ~1e17 raw, ~3e13
scaled), and no basis direction may be dropped — dropping breaks the symmetry of the Ritz
space and splits degenerate multiplets. A full run takes ~15 min.

**Validation** (`python modal_analysis.py --selftest`): polytope moments agree with exact
values (tesseract) and an independent cone quadrature (5-cell, 24-cell) to 1e-15; 4-volumes
match the literature (120-cell 787.857, 600-cell 26.4754, …). The solver reproduces the
analytic free-vibration frequencies of the unit 3-ball and 4-ball (torsional modes from
(l−1) J_ν(z) = z J_{ν+1}(z), ν = (n−2+2l)/2, and the radial breathing mode) to 1e-9…1e-15
(6e-6 for the higher l = 1 twist mode at degree 10), including the 4D multiplicities
2l(l+2) = 6, 16, 30. The lowest 4-ball mode is the 16-fold torsional l = 2 mode at
ωR/c_s = 2.6886.

**Accuracy.** Polytope edges are weak singularities of the elastic field, so convergence in
the degree is algebraic rather than exponential: the lowest multiplets are typically
accurate to 0.1–0.3 %. Each multiplet carries `rel_error` = relative change from degree 10
to 12; export stops at the first multiplet where this exceeds 1 % (`--max-error`), which
leaves 120–599 modes (26–106 multiplets) per polychoron, up to 2.5–3.6 × the lowest
frequency. Symmetry-degenerate modes agree to ~1e-9. Distinct multiplets closer than 1e-5
(relative) are merged into one: such near-degeneracies occur (3.7e-6 between a 9- and a
6-fold multiplet of prico), are inaudible, and their numerical eigenvectors mix, so only
the merged sums are exactly symmetric. Gains of symmetry-equivalent vertices and cells
agree to ≤ 5e-6.

### Output format `modal_output/<name>.json`

| field | meaning |
|---|---|
| `modes.frequency[g]` | dimensionless frequency of multiplet g (ascending); f_Hz = frequency · c_s/R |
| `modes.multiplicity[g]` | number of degenerate modes in the multiplet (symmetry) |
| `modes.rel_error[g]` | estimated relative frequency error |
| `gains.mean.normal[g]`, `.tangential[g]` | excitation of multiplet g by a hit at a uniformly random boundary point — the typical hit |
| `gains.cells[s]` | same, averaged over one cell of class s (`label` e.g. `truncated octahedron`, `cube #2`; `count` cells) |
| `gains.vertex` | same, for a hit exactly at a vertex (all vertices of a uniform polychoron are equivalent) |
| `cell_class[c]` | class s of cell c of `topology_output/<name>.json` (same indices) |
| `model` | physics, units, polynomial degree, error and merge tolerances |
| `geometry` | `volume` (in R⁴), `edge_length` (in R), element counts, detected symmetries |

A gain is mass · Σ_{k∈g} (u_k(p) · d)² for mass-normalised mode shapes u_k (∫ρ|u|² = 1): the
excitation of the multiplet by a unit impulse along d at p, picked up along d, relative to
the rigid-body value 1/mass (dimensionless, typically 0.01–10). `normal` uses d = outward
cell normal (vertex: normalised sum of the incident cell normals), `tangential` averages
over the 3 tangent directions. Summing over the multiplet makes a gain independent of the
arbitrary basis inside a degenerate eigenspace, hence identical for all symmetry-equivalent
cells — one table per cell class. Cell classes are orbits, not shapes: sidpith has two
classes of cubes (24 + 8). Sum and multiplicity are complementary: averaged over the whole
boundary (`gains.mean`), every mode of a symmetry multiplet is excited equally (Schur's
lemma), so the sum is `multiplicity` × the gain of one mode; at a single hit point the
modes share it unevenly and basis-dependently, and only the sum is physical. The ideal
object needs only the sum (all m modes ring as one sinusoid); the multiplicity matters once
the degeneracy is broken (below).

Why averages and not point values: at symmetric points (vertex, edge midpoint, face or
cell centre) many multiplets vanish by symmetry, so they are atypical hit points, and the
gain field of a multiplet peaks strongly at the cell's corners (up to ~150× its cell
average on the 5-cell). Interpolating point values therefore misses the spectrum of a
random hit by ~100 %; the cell averages are computed exactly (per-cell monomial moments
as quadratic forms, cross-checked against a degree-25 quadrature to 1e-11). In the 19
roundest polychora the lowest multiplet is a torsion mode (as in the 4-ball) whose
boundary motion is almost purely tangential: a perpendicular hit barely excites it
(< 5 % of the strongest multiplet), a glancing hit does; in the angular ones it is
excited normally as well (up to ~60 %: hex, pen).

**Synthesis recipe (game).** A hit with impulse J (normal component J_n, tangential J_t) on
cell c rings as a bank of damped oscillators, one per multiplet (s = `cell_class[c]`):

    y(t) = Σ_g A_g exp(-t / tau_g) cos(2π f_g t)
    f_g   = scale · frequency[g]                     (scale = c_s / R in Hz)
    A_g   = (J_n² cells[s].normal[g] + J_t² cells[s].tangential[g]) / mass · M(f_g)
    tau_g = 1 / (π · loss · f_g)                      (material loss factor)

M(f) is the spectrum of the contact force (mallet hardness; a half sine of contact time T
has |M(f)| = |cos(π f T)| / |1 − (2 f T)²|). Each oscillator is one two-pole resonator per
sample (a1 = 2 r cos(2π f/sr), a2 = −r², r = exp(−π · loss · f / sr)); hits add linearly.
For hits close to a vertex, blend towards `gains.vertex`; where only one sound per
polychoron is wanted, use `gains.mean`. If sound radiation is modelled separately, √gain
is the pure excitation amplitude.

**Symmetry breaking (detuning).** A real object is never perfectly symmetric: small
imperfections split every multiplet into m = `multiplicity` partials with slightly
different frequencies — the source of the beating ("warble") of bells. Random-matrix model
of an unknown imperfection: per multiplet, the relative detunings are the eigenvalues of a
random symmetric m×m matrix (GOE, which gives realistic level repulsion) scaled to the
chosen RMS, and the multiplet's gain is divided among the partials with Porter–Thomas
weights z_k²/|z|², z ~ N(0, I), so the attack stays the same. The split frequencies belong
to the object (fixed seed per exemplar), the weights are drawn per hit. RMS 0.02–0.2 % is
bell-like, ~0.3 % gives clearly audible slow beating, ~1 % chorus/roughness. The beating
pattern reflects the symmetry group: the 5-cell's multiplets have 1–6 modes, the H4
polychora's up to 48, which turn into dense shimmering clouds. In a game, capping each
multiplet at 3–4 partials sounds almost the same and keeps the oscillator count low.

### Example sounds

`python synth_modal.py` (or `make sounds`) renders one typical hit per polychoron to
`modal_output/sounds/<name>.wav` and all of them, highest to lowest, to `_tour.wav`
(legend with start times in `_tour.txt`). `--scale` (default 1000 Hz) is the common size
scaling c_s/R — e.g. steel (c_s ≈ 3200 m/s) with R = 3.2 m — and puts the fundamentals at
391 Hz (grip) … 521 Hz (hex). Other options: `--site mean|cell|vertex|all` (`all`: typical
hit, each cell class, vertex in sequence), `--direction normal|tangential`, `--loss`,
`--contact-ms` (default: mallet adapted to the exported band so its truncation is
inaudible), `--pickup displacement|velocity|acceleration`, `--duration`, `--tour-seconds`,
`--detune` (RMS detuning of the degenerate modes, default 0 = ideal object) with `--seed`
(which imperfect exemplar), e.g. `python synth_modal.py --detune 0.003 --out-dir
modal_output/sounds_detune`. The WAVs are git-ignored.

What sets the pitch of a solid polychoron: for one shape every size measure is equivalent
(f ∝ 1/size), across shapes the circumradius predicts it best — at equal R the lowest
frequencies of all 47 lie within ×1.33 (≈ 5 semitones), at equal 4-volume within ×2.2,
at equal edge length (as in `topology_output`) within ×22. The shape contributes a factor
between 0.46 (pentachoron: pointy, fills 3 % of its circumscribed ball, hence "soft") and
1.00 relative to the 4-ball of equal volume. The big H4 polychora (hi, ex, grix, grahi,
prahi, gidpixhi, …) ring within 0.5 % of that ball and share its partials (1 : 1.45 : 1.82 :
2.17 : 2.50 for multiplets 9, 16, 25, 36, …), so as solids they are practically
indistinguishable; the angular ones have their own partial patterns.

## Hollow bodies (thin-walled)

`python modal_hollow.py` (or `make hollow`) computes the same kind of modal data for the
polychora as **thin-walled hollow 4D bodies** → `modal_hollow_output/<name>.json`, same
format as `modal_output` (differences below). As solids, polychora of equal circumradius
sound alike; as hollow bodies they sound like their walls, whose pitch is set by the cell
shapes and sizes — the lowest frequencies of the 47 spread over ×16 (0.124 for grip to 1.985
for ex at h_ref) instead of ×1.33, and e.g. the 120-cell and the 600-cell, indistinguishable
as solids, are 1.7 octaves apart.

**Model.** Every cell is a flat 3D wall ("hyperplate") of thickness h that bends in the 4th
direction, along its normal (Kirchhoff: (h³/12) ∫ [λ*(Δw)² + 2μ|∇∇w|²], λ* = 2λμ/(λ+2μ),
plane stress σ_nn = 0). Thin-wall limit: a ridge (polygon shared by two cells) cannot move,
because any motion in its 2D normal plane stretches at least one wall in-plane, which is
infinitely stiffer than bending as h → 0; so each wall's deflection vanishes on its faces.
Neighbouring walls are welded: the hinge rotation across a ridge (the slope, taken in the
ridge's normal plane with the orientation sign det[[n_A·n_B, n_A·m_B], [m_A·n_B, m_A·m_B]])
is continuous, imposed by a penalty β = 10³ μh³/ℓ. Same material and units as the solids,
reference thickness h_ref = 0.06 R; the bending spectrum is exactly linear in h:
**f_Hz = frequency · (h/h_ref) · c_s/R**.

**Method.** Per cell, w = b(s)·poly(s) with a "flat-top" bubble b = Π_faces tanh(dist_f/δ):
it vanishes linearly on every face like the polynomial bubble Π dist_f (its δ → ∞ limit),
but stays ≈ 1 inside even for the 62-face cells, where the polynomial bubble collapses to a
narrow peak. The p → ∞ limit does not depend on δ, the speed does: δ = 3 × inradius for cells
with ≤ 8 faces (nearly the polynomial bubble — tetrahedral walls converge to 0.1 % at degree
6 instead of 1.4 % with δ = r/2), 1 × for ≤ 12 faces, 0.5 × for more. Cell matrices are built
once per cell shape (Gauss quadrature on the flag tetrahedra) and carried to every
congruent cell by an orthogonal map found from the vertex sets. Degrees are chosen per
polychoron within a budget of 26 000 unknowns: the soft (large)
walls get degree 4–6, much stiffer small walls 1–2 (in the audible band they only transmit
rotations). Modes are computed up to 3 × the lowest frequency (at least 6 per soft wall,
≤ 800), and the export ends with a complete band: the computed set or the 3× limit may end
inside a band, whose modes would then be under-represented, so the last band is cut at a
band gap. This leaves 10–720 modes (2–59 multiplets) per polychoron, up to 1.5–2.7 × the
lowest frequency (the isolated-wall polychora below: 6–9 bands up to 4.0–5.4 ×). A full run
takes ~2.5 h.

**Isolated walls.** In prahi, prix and gidpixhi the 120 largest walls touch no other large
wall, only smaller, stiffer ones: the lowest mode any neighbour wall can have (simply
supported) lies at 4.6–5.6 × the fundamental. Every band is then one mode of a single large
wall, repeated on all 120 of them: one multiplet per band, multiplicity 120 × wall multiplet.
`modal_hollow.py` detects this (a single soft shape, no ridge between two soft walls,
neighbours ≥ `--isolated-min` 3.5 × f1) and computes one large wall welded to all its 32–62
neighbours, with the smooth-distance basis below (degree 10, neighbours 8). With the
neighbours' outer ridges clamped, all other large walls are at rest, which to first order is
the band's mean (its trace). Checked against the full coupled prix at low degree (50 280
unknowns, 22 min, 5 GB): band 1 spans 0.870–0.941 with mean 0.936, the patch gives 0.941.
Hinging the outer ridges towards the other large walls moves into the lower part of the band;
that difference (≤ 1.2 %) is a rough estimate of the band's half width and enters `rel_error`,
together with the changes from wall degree 8 → 10 and neighbour degree 6 → 8. The export ends
at the neighbours' lowest simply supported mode, below which no band of theirs can lie (above
it, their bands of 600–1200 modes would dominate). Result: 6–9 bands up to 4.0–5.4 × f1
instead of 2 up to 1.9 ×, `rel_error` ≤ 1.9 %. The welds are not rigid: compared with a
clamped single wall, the neighbours lower the frequencies by about 4 % (prahi), 7 %
(gidpixhi) and 11 % (prix). The neighbours' `gains.cells` hold the small share of each band
that leaks into them (walls that touch no large wall: 0). grix, rox, srix, sidpixhi and prit
also have isolated largest walls, but their neighbours resonate already at 1.7–2.1 × f1; they
stay with the coupled model (`--no-isolated` forces it everywhere).

**Validation and accuracy.** A single simply supported cube wall reproduces the analytic
f = 3π² h √((λ*+2μ)/12) / (2π) to 5 digits (and its (1,1,2) multiplet at exactly 2×);
frequencies scale exactly with h; all multiplets follow the symmetry groups; on the
tesseract, degree 9 with δ = r/3, r, 3r gives the same f1 = 0.1674 and multiplets, in
agreement with an independent polynomial-bubble prototype. A first approach with free
polynomial displacement fields per cell and displacement penalties locked (frequencies
independent of h) — the fixed-ridge basis avoids that. Walls with obtuse inner dihedral
angles (dodecahedron 117°, Archimedean cells up to ~160°) have weak edge singularities, so
convergence in the degree is algebraic: first bands are accurate to ~1 % for cube and
tetrahedron walls and ~5–10 % (too high; Ritz values are upper bounds) for dodecahedral and
many-faced walls — probably more, see the next paragraph. `rel_error` is an indicator per
multiplet (the clamped single wall of the dominant cell shapes, degree p vs p + 2),
pessimistic for higher modes; for the first multiplet it is 2.5 % (median) and at most 7 %.

**Smooth-distance basis.** The flat-top bubble converges slowly when a wall is clamped or
stiffly welded: its slope at a face is poly/δ there, so the polynomial itself would have to
vanish on all 26–62 faces. The isolated-wall model uses φ = (Σ_f d_f⁻²)^(−1/2) instead, a
smooth distance to the face planes (φ ≈ d_f at each face, no layer width): φ²·poly has zero
value and slope on every face and carries the interior, φ·poly (degree p − 2) the rotations
at the welds; each cell's basis is orthonormalised. A clamped rhombicuboctahedron wall
converges to 0.03 % at degree 4 (flat-top bubble: 2 % at degree 10), the large walls of
prahi, prix and gidpixhi to 0.1–0.5 % at degree 8 for their lowest four modes. Measured
against it, the flat-top results for these three were too high by 33 % (prahi), 16 % (prix)
and 31 % (gidpixhi), far beyond their indicator (7, 4, 6 %): the change p → p + 2 of a
slowly converging sequence underestimates the distance to its limit, and the stiff
neighbours at degree 1–2 could hardly rotate. The other polychora still use the flat-top
bubble; walls with many faces are probably too high there as well, ratios within one
polychoron less so.

**Differences to `modal_output`.** No `gains.vertex` (vertices lie on the fixed ridges),
`tangential` = 0 (in-plane membrane motion is not modelled; its modes are far higher),
`gains.mean` = multiplicity (every mode bends all the mass, so a random hit excites each
mode equally on average), and `gains.cells` depend strongly on the hit wall: a hit excites
mainly the bands of its own wall shape, a small stiff wall barely excites the low bands.
`cell_class` is the same as in `modal_output`.

**Character and limits.** A hollow polychoron's spectrum consists of bands: one mode per
soft wall and wall mode, split by the coupling through the welds. Where the largest walls
do not touch each other (prahi, prix, gidpixhi, …) the coupling runs only through smaller,
stiffer walls and a band of 120 modes stays narrow — the body sounds like a single wall. The
three isolated-wall polychora therefore share their first four tones (0, ≈ 1025, ≈ 1773,
≈ 2042 cent: the clamped modes of a nearly spherical wall) and differ in pitch (f1 = 0.995,
0.837, 0.805) and in their upper tones.
`synth_modal.py --detune` makes the clusters shimmer, as real imperfections would. Only the lowest bands are computed,
so the high "ping" of a small stiff wall hit directly is missing. Example sounds:
`make hollow-sounds` (`synth_modal.py --modal-dir modal_hollow_output --out-dir
modal_hollow_output/sounds`; `--scale` = c_s/R for the reference thickness).

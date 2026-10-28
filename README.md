# arostica-benchmark-cm

Verification of the dynamic solver of [CardioMechanics](https://github.com/tobias-gerach/CardioMechanics/tree/experimental)
on the cardiac elastodynamics benchmark of Aróstica et al. (2025),
[doi:10.1016/j.cma.2024.117485](https://doi.org/10.1016/j.cma.2024.117485): the displacement at
probe points over one heart beat, scored by RED, eq. (21), against the cohort of codes that took part.
Every run uses CardioMechanics' quadratic tetrahedron, T10. The commits the CSVs record are on the
`experimental` branch of that fork.

The report of the final runs is published on GitHub Pages from `docs/index.html`: the repository's
Pages source is the `main` branch, folder `/docs`. `html_report.py` rebuilds it from the benchmark
data, which CI does not have, so the built file is committed.

| File | |
|---|---|
| `mono.py` | Benchmark 1, the monoventricle: fibers, mesh, runs, scoring |
| `biv.py` | Benchmark 2, the biventricle on the coarse and fine mesh: meshes, runs, scoring |
| `common.py` | conversion to tetgen, fiber bases, deck, runner, VTK probing, cohort, RED, plots, records |
| `dolfin_fibers.py` | runs in the toolkit's dolfin image and writes the fibers at the degree-2 nodes |
| `templates/mono.xml` | the deck of Tables 1 to 3 |
| `templates/biv.xml` | the deck of Tables 2, 7, 8 and 9 |
| `html_report.py`, `templates/report.html` | the self-contained HTML report `docs/index.html` of the final runs, served by GitHub Pages |
| `meshes/mono/` | the converted T10 mesh and its bases, committed |
| `meshes/biv_coarse/`, `meshes/biv_fine/` | the converted T10 meshes and their bases, rebuilt by `biv.py mesh` and not committed |
| `results/` | one CSV row per case in `mono.csv` and per mesh in `biv.csv`, the probe traces, the overlay plots and the pass/fail tables `mono_report.md` and `biv_report.md` |

## Usage

`DATA` is the unpacked benchmark dataset, `BIN` the directory of the CardioMechanics binaries.

```sh
python -m venv .venv && .venv/bin/pip install -r requirements.txt
# once: the fibers by the toolkit's own Docker image, then the mesh and bases (committed)
.venv/bin/python mono.py fibers --data $DATA --out $OUT
.venv/bin/python mono.py mesh --data $DATA --bin $BIN --fibers $OUT/fibers_p2.npz
# a run; raw output goes to $OUT/runs and is not committed
.venv/bin/python mono.py run --data $DATA --bin $BIN --out $OUT --case 0A --np 16
# the pass/fail table of every case in results/mono.csv
.venv/bin/python mono.py report
```

The biventricle takes the meshes and fibers of the data as they are:

```sh
.venv/bin/python biv.py mesh --data $DATA --bin $BIN
.venv/bin/python biv.py run --data $DATA --bin $BIN --out $OUT --mesh coarse --np 16
.venv/bin/python biv.py report
# the HTML report of both benchmarks, docs/index.html
.venv/bin/python html_report.py --data $DATA
```

Cases are 0A, 0B, 1, A, B and C. The biventricle's meshes are coarse (Ω_h1) and fine (Ω_h2).

The branch and commit of a run are those of the git checkout holding `BIN`.

## Method

- **Mesh**: the distributed `ellipsoid_0.005.h5` is written as a tetgen T4 mesh with surfaces
  endocardium 1, epicardium 2 and base 3, normals pointing out of the wall, and converted to T10 by
  CardioMechanics' `ConvertT4toT10`. The T4 mesh is an intermediate and is not kept. The mid-edge
  nodes sit at the edge midpoints.
- **Fibers**: the toolkit solves the transmural Laplace problem and interpolates f, s, n of eq. (15)
  at degree 2. Each T10 node takes the values at its point, and the T10 shape functions interpolate
  them to the element's centroid, which holds the element's single frame. f is
  normalized, s orthogonalized against f, and n = f × s. The paper's n, eq. (15), is s × f, the
  opposite sign; the Holzapfel law does not use n. The toolkit's n is transmural and its s lies in
  the wall, across the fibers, unlike CardioMechanics' usual transmural sheet. The Holzapfel law
  takes s from the second vector of each basis, so it sees the toolkit's s, as the cohort's codes
  and Table 2 do.
- **Deck**: Holzapfel with k = 100 and κ = 1e6, Kelvin-Voigt η = 100 Pa·s, ρ = 1000 kg/m³, Bestel
  with PK2 tension, Robin epicardium normal-only (α 1e8, β 5e3) and base (α 1e5, β 5e3),
  GeneralizedAlpha at RhoInf = 0.5, dt = 1 ms over 1 s, full Newton with a direct MUMPS solve,
  quadrature degree 5, VTU every 10 ms. The endocardial pressure of Table 4 is a
  `BestelPressure` of `ApplyPressureFromFunction`. Step 0A switches it off, Step 0B sets the
  contractility to zero, Step 1 keeps both, and Cases A to C take the stiffness and contractility of
  Table 5.
- **Pressure sign**: eq. (1) writes the endocardial traction as +p J F^-T N. The cohort's Step 0B
  inflates the cavity for positive p: at p0, on the midwall at y = 0.03, every code moves out
  (uy > 0) and, at the peak, towards the apex (ux < 0). So the paper's N on the endocardium points
  out of the cavity, into the wall. `ApplyPressureFromFunction` applies the traction p n along the
  face normal, which for these meshes points out of the wall, into the cavity, so the deck scales
  `BestelPressure` by -1.
- **Probing**: the displacement at p0 = (0.025, 0.03, 0) and p1 = (0, 0.03, 0) is interpolated with
  the shape functions of the VTK cell holding the point in the reference configuration.
- **RED**: eq. (21) against the index-wise mean of the ten codes of the paper's figures, matching
  samples by index as the paper does. This reproduces Table 11. A case passes if CardioMechanics' RED
  at every probe point is no worse than the worst cohort member's, computed the same way. For Cases A
  to C these thresholds reproduce the worst RED of Table 11, which `mono.py report` asserts. A run
  that fails is recorded with the time it failed at and fails its case.

## Method: biventricle

- **Mesh**: `bi_ventricular.h5` of `biv_mesh_and_fibers_{coarse,fine}_p2` is converted as the
  monoventricle's, keeping the data's tags as surfaces: base 10, right endocardium 20, left
  endocardium 30 and epicardium 40. The left endocardium is the round cavity around the x axis, the
  right the crescent at z > 0, as eqs. (19) and (20) place them.
- **Fibers**: `fibers_biv.vtk` holds LifeX's f0, s0 and n0 at the P2 nodes, the vertices and edge
  midpoints. Each T10 node takes the values at its point, matched by coordinate: the file rounds the
  points to six digits, up to 1e-6 m off, against 4e-4 m between the nearest nodes, and every node
  must find a point of its own. The frames are then built as the monoventricle's, with s0, the
  toolkit's `sheet`, as s: the data already carries the paper's convention, which footnote 6 says
  swaps s and n against the rule-based method's. LifeX's n0 equals f0 × s0.
- **Deck**: Table 2's Holzapfel parameters and Table 7: κ = 1e6, k = 100, η = 100 Pa·s,
  ρ = 1000 kg/m³, Robin epicardium normal-only (α 1e8, β 5e3) and base (α 1e6, β 5e3). Bestel as in
  Table 8 with σ0 = 1.5e5, t_sys = 0.163 s and t_dias = 0.5 s. Each endocardium carries its own
  `BestelPressure` group of Table 9, scaled by -1 as the monoventricle's. RhoInf 0.5, quadrature
  degree 5, the 1 ms step and the Newton and MUMPS settings are the monoventricle's.
- **Probing and RED**: as for the monoventricle, with p2 = (0.025, 0, 0.072) added, against the
  index-wise mean of the nine codes of Table 12. The `simvascular` pickle is SimVascular's P1 result.
- **Table 12**: the coarse thresholds computed from the pickles reproduce it, and `biv.py report`
  asserts so. The fine ones do not, because the fine CARPentry pickle is a copy of the coarse one,
  every sample equal; the paper used CARPentry's fine result. The copy enters the cohort mean as well
  as the threshold, so the fine RED is scored against a mean that differs slightly from the paper's.
  The verdict uses the pickles, which give thresholds of 0.997, 0.562 and 1.891 at p0, p1 and p2
  against Table 12's 1.019, 0.513 and 1.769.
- **Thresholds**: CARPentry at p0 and Simula at p1 and p2 lie far from the other codes on both
  meshes, with RED up to 1.9. SimVascular P1 follows, up to 0.95 at p2, and the six other codes stay
  below 0.42. A pass against the worst code is therefore weak evidence, and the six close codes are
  the more telling comparison.

## Data and license

The benchmark dataset of Aróstica et al. is
[doi:10.5281/zenodo.14260459](https://doi.org/10.5281/zenodo.14260459), under CC BY 4.0. This
repository redistributes material derived from it: the monoventricle mesh and fibers in
`meshes/mono/`, made with the dataset's toolkit (MIT License, copyright (c) 2023 Reidmen Arostica),
and the cohort's displacement curves embedded in `docs/index.html`.

The code of this repository is under the MIT License, see `LICENSE`.

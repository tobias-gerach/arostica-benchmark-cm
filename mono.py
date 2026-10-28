#!/usr/bin/env python3
"""Benchmark 1 of Aróstica et al. (2025): the monoventricle.

    mono.py fibers --data DATA --out OUT          degree-2 fibers by the toolkit's Docker image
    mono.py mesh --data DATA --bin BIN --fibers NPZ   T10 mesh and bases into meshes/mono
    mono.py run --data DATA --bin BIN --out OUT --case 0A --np 16
    mono.py report                                 the pass/fail table of results/mono.csv

DATA is the benchmark data repository, BIN the directory of the CardioMechanics binaries.
"""
import argparse
import csv
import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np

import common

MESH = common.REPO / "meshes" / "mono"
STEM = "ellipsoid"
H5 = Path("meshes_and_fibers/data/monoventricle_mesh/ellipsoid_0.005.h5")
IMAGE = "cardiac-benchmark-toolkit"
PROBES = {"p0": (0.025, 0.03, 0.0), "p1": (0.0, 0.03, 0.0)}
HOLZAPFEL = {"a": 59.0, "af": 18472.0, "afs": 216.0, "as": 2481.0}       # Table 2
# Per case: the pickle prefix of the cohort and the parameters that differ between cases. pressure
# switches the endocardial pressure of Table 4 on (1) or off (0). Cases A to C take the stiffness and
# contractility of Table 5.
CASES = {"0A": ("monoventricular_nonblinded_step_0A", {**HOLZAPFEL, "sigma0": 1.5e5, "pressure": 0}),
         "0B": ("monoventricular_nonblinded_step_0B", {**HOLZAPFEL, "sigma0": 0, "pressure": 1}),
         "1": ("monoventricular_nonblinded_step_1", {**HOLZAPFEL, "sigma0": 1.5e5, "pressure": 1}),
         "A": ("monoventricular_blinded_A", {"a": 177.0, "af": 55416.0, "afs": 648.0, "as": 7443.0, "sigma0": 2e5, "pressure": 1}),
         "B": ("monoventricular_blinded_B", {"a": 295.0, "af": 92360.0, "afs": 1080.0, "as": 12405.0, "sigma0": 1e5, "pressure": 1}),
         "C": ("monoventricular_blinded_C", {"a": 19.0, "af": 6157.0, "afs": 72.0, "as": 827.0, "sigma0": 2e5, "pressure": 1})}
# The worst cohort RED of Table 11 at p0 and p1, which the thresholds computed from the pickles must
# reproduce. Ambit's 0.79 at p0 of Case B is 0.179 in the pickles and not the worst.
TABLE11 = {"A": (0.276, 0.371), "B": (0.360, 0.352), "C": (0.156, 0.211)}
RHO_INF = 0.5                                                           # see templates/mono.xml
QUADRATURE = 5                                                          # see templates/mono.xml


def fibers(args):
    """Build the toolkit's image from the data repository and run dolfin_fibers.py in it."""
    toolkit = args.data / "meshes_and_fibers"
    subprocess.run(["docker", "build", "-q", "-t", IMAGE, "-f", "docker/Dockerfile", "."], cwd=toolkit, check=True)
    args.out.mkdir(parents=True, exist_ok=True)
    # dolfin compiles its forms into HOME, which the image does not provide. The image's OpenBLAS
    # predates Sapphire Rapids, and the kernels it picks there give NaN from the toolkit's LU solve;
    # the Skylake-X kernels use the same AVX-512 instructions.
    subprocess.run(["docker", "run", "--rm", "--user", str(os.getuid()), "--tmpfs", "/tmp:exec,mode=1777", "-e", "HOME=/tmp",
                    "-e", "OPENBLAS_CORETYPE=SkylakeX",
                    "-v", f"{args.data.resolve()}:/data:ro", "-v", f"{common.REPO}:/repo:ro", "-v", f"{args.out.resolve()}:/out",
                    IMAGE, "python3", "/repo/dolfin_fibers.py", f"/data/{H5.with_suffix('.xdmf')}", "/out/fibers_p2.npz"],
                   check=True)


def mesh(args):
    MESH.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as t4:
        common.write_tetgen_t4(args.data / H5, Path(t4), STEM)
        common.convert_t10(args.bin / "ConvertT4toT10", Path(t4) / STEM, MESH / STEM)
    with np.load(args.fibers) as fibers:
        common.write_bases(fibers["X"], fibers["f"], fibers["s"], MESH, STEM)


def run(args):
    """Run a case and record it. A run that fails is recorded with the time it failed at and no RED."""
    prefix, values = CASES[args.case]
    wd = args.out / "runs" / "mono" / args.case
    (wd / "out").mkdir(parents=True, exist_ok=True)
    common.write_deck("mono.xml", wd, step=args.case, mesh=MESH, rho_inf=RHO_INF, quadrature=QUADRATURE, **values)
    print(f"mono {args.case} ...", flush=True)
    wall, finished = common.run(args.bin / "CardioMechanics", wd, args.np)
    print(f"{wall:.0f} s", flush=True)

    curves = common.cohort(args.data / "results_time_curves" / "data", prefix, PROBES)
    outcome = common.evaluate(wd, finished, MESH / f"{STEM}.node", PROBES, curves, f"mono_{args.case}",
                              f"Monoventricle Step {args.case}, T10, RhoInf {RHO_INF:g}")
    common.record(common.RESULTS / "mono.csv", {
        "case": args.case, "dofs": common.dofs(MESH / f"{STEM}.ele"), "np": args.np, "rho_inf": f"{RHO_INF:g}",
        "quadrature": QUADRATURE, **outcome, "newton": common.newton_iterations(wd), "wall_s": f"{wall:.0f}",
        **common.provenance(args.bin)})


def report(args):
    """Write results/mono_report.md: the latest run of every case and its verdict per probe point."""
    with open(common.RESULTS / "mono.csv") as f:
        runs = {r["case"]: r for r in csv.DictReader(f)}
    header = ["Case", "DOFs", *(f"RED {p}" for p in PROBES), *(f"cohort max {p}" for p in PROBES),
              *(f"Table 11 {p}" for p in PROBES), "Verdict", "Commit", "Plot"]
    lines = ["# Monoventricle: pass/fail", "",
             f"T10, RhoInf {RHO_INF:g}, quadrature degree {QUADRATURE}. A run passes if its RED at every probe point is no "
             "worse than the worst cohort member's (cohort max). The verdict names the probe points that fail. Cases A to C "
             "list the worst RED of Table 11 next to it.", "",
             "| " + " | ".join(header) + " |", "|---" * len(header) + "|"]
    for case, run in sorted(runs.items(), key=lambda item: list(CASES).index(item[0])):
        failed = run.get("failed_at")
        thresholds = [float(run[f"cohort_max_red_{p}"]) for p in PROBES]
        if case in TABLE11:
            assert np.allclose(thresholds, TABLE11[case], rtol=0, atol=5e-4), f"{case}: {thresholds} != Table 11"
            table11 = [f"{t:.3f}" for t in TABLE11[case]]
        else:
            table11 = ["–"] * len(PROBES)
        if failed:
            verdict = f"fail at t = {failed} s"
        else:
            worse = [p for p, t in zip(PROBES, thresholds) if float(run[f"red_{p}"]) > t]
            verdict = f"fail ({', '.join(worse)})" if worse else "pass"
        cells = [case, run["dofs"], *(run.get(f"red_{p}") or "–" for p in PROBES), *(f"{t:.4f}" for t in thresholds),
                 *table11, verdict, run["commit"], "–" if failed else f"[png](mono_{case}.png)"]
        lines.append("| " + " | ".join(cells) + " |")
    (common.RESULTS / "mono_report.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("report")
    for name in ("fibers", "mesh", "run"):
        command = commands.add_parser(name)
        command.add_argument("--data", type=Path, required=True, help="the benchmark data repository")
    commands.choices["fibers"].add_argument("--out", type=Path, required=True)
    for name in ("mesh", "run"):
        commands.choices[name].add_argument("--bin", type=Path, required=True, help="directory of the CardioMechanics binaries")
    commands.choices["mesh"].add_argument("--fibers", type=Path, required=True, help="fibers_p2.npz of the fibers command")
    commands.choices["run"].add_argument("--out", type=Path, required=True, help="directory of the raw output")
    commands.choices["run"].add_argument("--case", choices=CASES, default="0A")
    commands.choices["run"].add_argument("--np", type=int, default=16)
    args = parser.parse_args()
    {"fibers": fibers, "mesh": mesh, "run": run, "report": report}[args.command](args)

#!/usr/bin/env python3
"""Benchmark 2 of Aróstica et al. (2025): the biventricle, on the coarse and the fine mesh.

    biv.py mesh --data DATA --bin BIN              T10 meshes and bases into meshes/biv_*
    biv.py run --data DATA --bin BIN --out OUT --mesh coarse --np 16
    biv.py report                                  the pass/fail table of results/biv.csv

DATA is the benchmark data repository, BIN the directory of the CardioMechanics binaries. The meshes
are rebuilt from DATA in seconds and are not committed.
"""
import argparse
import csv
import tempfile
from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

import common
import mono

STEM = "biv"
BASE, ENDO_RV, ENDO_LV, EPI = 10, 20, 30, 40                    # surface tags of the biventricle
PROBES = {**mono.PROBES, "p2": (0.025, 0.0, 0.072)}
# The worst cohort RED of Table 12 at p0, p1 and p2. The thresholds computed from the pickles
# reproduce the coarse row; the fine pickle of CARPentry is a copy of its coarse one, see the README.
TABLE12 = {"coarse": (0.915, 0.848, 1.472), "fine": (1.019, 0.513, 1.769)}


def data_dir(data, resolution):
    return data / "meshes_and_fibers" / "data" / f"biv_mesh_and_fibers_{resolution}_p2"


def mesh_dir(resolution):
    """The T10 mesh of resolution and its bases."""
    return common.REPO / "meshes" / f"biv_{resolution}"


def mesh(args):
    """Convert both meshes and give each its LifeX P2 fibers. The provided sheet is s: the paper's
    footnote 6 swaps s and n against the rule-based method's own convention."""
    for resolution in TABLE12:
        source, target = data_dir(args.data, resolution), mesh_dir(resolution)
        target.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory() as t4:
            common.write_tetgen_t4(source / "bi_ventricular.h5", Path(t4), STEM, (BASE, ENDO_RV, ENDO_LV, EPI))
            common.convert_t10(args.bin / "ConvertT4toT10", Path(t4) / STEM, target / STEM)
        reader = vtk.vtkGenericDataObjectReader()
        reader.SetFileName(str(source / "fibers_biv.vtk"))
        reader.ReadAllVectorsOn()
        reader.Update()
        grid = reader.GetOutput()
        f, s = (vtk_to_numpy(grid.GetPointData().GetArray(k)).astype(float) for k in ("f0", "s0"))
        common.write_bases(vtk_to_numpy(grid.GetPoints().GetData()).astype(float), f, s, target, STEM)


def run(args):
    """Run a mesh and record it. A run that fails is recorded with the time it failed at and no RED."""
    wd = args.out / "runs" / "biv" / args.mesh
    mesh_path = mesh_dir(args.mesh)
    (wd / "out").mkdir(parents=True, exist_ok=True)
    common.write_deck("biv.xml", wd, mesh=mesh_path, resolution=args.mesh, rho_inf=mono.RHO_INF,
                      quadrature=mono.QUADRATURE, **mono.HOLZAPFEL)
    print(f"biv {args.mesh} ...", flush=True)
    wall, finished = common.run(args.bin / "CardioMechanics", wd, args.np)
    print(f"{wall:.0f} s", flush=True)

    curves = common.cohort(args.data / "results_time_curves" / "data", f"biventricular_{args.mesh}", PROBES, common.BIV_COHORT)
    outcome = common.evaluate(wd, finished, mesh_path / f"{STEM}.node", PROBES, curves, f"biv_{args.mesh}",
                              f"Biventricle, {args.mesh} mesh, T10, RhoInf {mono.RHO_INF:g}")
    common.record(common.RESULTS / "biv.csv", {
        "mesh": args.mesh, "dofs": common.dofs(mesh_path / f"{STEM}.ele"), "np": args.np, "rho_inf": f"{mono.RHO_INF:g}",
        "quadrature": mono.QUADRATURE, **outcome, "newton": common.newton_iterations(wd), "wall_s": f"{wall:.0f}",
        **common.provenance(args.bin)})


def report(args):
    """Write results/biv_report.md: the latest run of every mesh and its verdict per probe point."""
    with open(common.RESULTS / "biv.csv") as f:
        runs = {r["mesh"]: r for r in csv.DictReader(f)}
    header = ["Mesh", "DOFs", *(f"RED {p}" for p in PROBES), *(f"cohort max {p}" for p in PROBES),
              *(f"Table 12 {p}" for p in PROBES), "Verdict", "Commit", "Plot"]
    lines = ["# Biventricle: pass/fail", "",
             f"T10, RhoInf {mono.RHO_INF:g}, quadrature degree {mono.QUADRATURE}. A run passes if its RED at every probe point "
             "is no worse than the worst cohort member's (cohort max), computed from the pickles. The verdict names the "
             "probe points that fail. The worst RED of Table 12 is listed next to it.", "",
             "| " + " | ".join(header) + " |", "|---" * len(header) + "|"]
    for resolution, run in sorted(runs.items(), key=lambda item: list(TABLE12).index(item[0])):
        failed = run.get("failed_at")
        thresholds = [float(run[f"cohort_max_red_{p}"]) for p in PROBES]
        if resolution == "coarse":
            assert np.allclose(thresholds, TABLE12[resolution], rtol=0, atol=5e-4), f"{thresholds} != Table 12"
        if failed:
            verdict = f"fail at t = {failed} s"
        else:
            worse = [p for p, t in zip(PROBES, thresholds) if float(run[f"red_{p}"]) > t]
            verdict = f"fail ({', '.join(worse)})" if worse else "pass"
        cells = [resolution, run["dofs"], *(run.get(f"red_{p}") or "-" for p in PROBES),
                 *(f"{t:.4f}" for t in thresholds), *(f"{t:.3f}" for t in TABLE12[resolution]), verdict, run["commit"],
                 "-" if failed else f"[png](biv_{resolution}.png)"]
        lines.append("| " + " | ".join(cells) + " |")
    (common.RESULTS / "biv_report.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("report")
    for name in ("mesh", "run"):
        command = commands.add_parser(name)
        command.add_argument("--data", type=Path, required=True, help="the benchmark data repository")
        command.add_argument("--bin", type=Path, required=True, help="directory of the CardioMechanics binaries")
    commands.choices["run"].add_argument("--out", type=Path, required=True, help="directory of the raw output")
    commands.choices["run"].add_argument("--mesh", choices=TABLE12, default="coarse")
    commands.choices["run"].add_argument("--np", type=int, default=16)
    args = parser.parse_args()
    {"mesh": mesh, "run": run, "report": report}[args.command](args)

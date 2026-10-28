#!/usr/bin/env python3
"""A self-contained HTML report of the final T10 runs of both benchmarks.

    html_report.py --data DATA        writes docs/index.html, which GitHub Pages serves

DATA is the benchmark data repository, which holds the cohort's curves. The runs are the latest rows
of results/mono.csv and results/biv.csv, with their probe traces.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

import biv
import common
import mono

TITLE = "CardioMechanics T10 on the Aróstica et al. (2025) cardiac elastodynamics benchmark"
# Maxima of the BestelPressure curves in Pa: the monoventricle's endocardium (Table 4), and the
# biventricle's left and right endocardium (Table 9), from the converged solution of ODE (7).
PEAK_MONO, PEAK_LV, PEAK_RV = 16117.52, 16491.15, 4171.07
MMHG = 133.322


def pressure(peak):
    return f"{peak / 1e3:.1f} kPa ({peak / MMHG:.0f} mmHg)"


def describe_mono(case, values):
    """What a monoventricle case loads the ventricle with, and how its stiffness and contractility
    differ from Step 1's."""
    tension = f"active tension with σ0 = {values['sigma0'] / 1e3:g} kPa" if values["sigma0"] else "no active tension"
    load = f"endocardial pressure up to {pressure(PEAK_MONO)}" if values["pressure"] else "no cavity pressure"
    # Table 5 scales a, af, afs and as by one factor; a is rounded to whole Pa, so af carries it best.
    factor = values["af"] / mono.HOLZAPFEL["af"]
    stiffness = "the passive stiffness of Table 2" if factor == 1 else f"{factor:.2g} times the passive stiffness of Table 2"
    kind = {"0A": "Active response only", "0B": "Passive response only", "1": "Active and passive response"}.get(
        case, "Blinded variation of Step 1")
    return f"{kind}: {tension}, {load}, {stiffness}."


def describe_biv(resolution, elements):
    return (f"Blinded biventricle on the {resolution} mesh ({elements:,} tetrahedra): active tension with σ0 = 150 kPa "
            f"(onset of systole 0.163 s, of diastole 0.5 s), left endocardial pressure up to {pressure(PEAK_LV)} and "
            f"right up to {pressure(PEAK_RV)}, the passive stiffness of Table 2.")


def latest(path, key):
    """The last row of path for each value of its column key."""
    with open(path) as f:
        return {r[key]: r for r in csv.DictReader(f)}


def entry(data, benchmark, label, description, name, prefix, probes, codes, row):
    """One run: its RED and rank among the cohort per probe point, the cohort's RED, and the curves in mm.
    A failed run has the time it failed at, and neither RED nor curve."""
    curves = common.cohort(data, prefix, probes, codes)
    mean = np.mean(list(curves.values()), axis=0)
    failed_at = float(row["failed_at"]) if row.get("failed_at") else None
    red = u_mm = None
    if failed_at is None:
        u = np.loadtxt(common.RESULTS / f"{name}_trace.csv", delimiter=",", skiprows=1)[:, 1:].reshape(common.SAMPLES, len(probes), 3)
        score = common.red(u, mean)
        # The traces are the ones the CSV row was scored from.
        assert np.allclose(score, [float(row[f"red_{p}"]) for p in probes], rtol=0, atol=5e-5), f"{name}: trace and CSV disagree"
        red, u_mm = score.round(4).tolist(), (u * 1e3).round(5).transpose(1, 2, 0).tolist()   # [probe][component][sample]
    cohort_red = {code: common.red(v, mean) for code, v in curves.items()}
    return {
        "benchmark": benchmark, "label": label, "description": description,
        "probes": list(probes), "red": red, "failed_at": failed_at,
        "cohort": {code: r.round(4).tolist() for code, r in cohort_red.items()},
        "dofs": int(row["dofs"]), "newton": int(row["newton"]), "wall_s": int(row["wall_s"]), "commit": row["commit"],
        "passed": row["pass"] == "True", "u": u_mm,
        "codes": {code: (v * 1e3).round(5).transpose(1, 2, 0).tolist() for code, v in curves.items()},
    }


def main(data):
    runs = []
    mono_rows = latest(common.RESULTS / "mono.csv", "case")
    biv_rows = latest(common.RESULTS / "biv.csv", "mesh")
    for case, (prefix, values) in mono.CASES.items():
        label = f"Step {case}" if case[0].isdigit() else f"Case {case}"
        runs.append(entry(data, "Monoventricle", label, describe_mono(case, values), f"mono_{case}", prefix,
                          mono.PROBES, common.COHORT, mono_rows[case]))
    for resolution in biv.TABLE12:
        elements = sum(1 for _ in open(biv.mesh_dir(resolution) / f"{biv.STEM}.ele")) - 1
        runs.append(entry(data, "Biventricle", f"{resolution.capitalize()} mesh", describe_biv(resolution, elements),
                          f"biv_{resolution}", f"biventricular_{resolution}", biv.PROBES, common.BIV_COHORT,
                          biv_rows[resolution]))
    meta = {"title": TITLE, "rho_inf": mono.RHO_INF, "quadrature": mono.QUADRATURE}
    html = (common.REPO / "templates" / "report.html").read_text()
    html = html.replace("__DATA__", json.dumps({"meta": meta, "runs": runs}, separators=(",", ":")))
    (common.REPO / "docs").mkdir(exist_ok=True)
    (common.REPO / "docs" / "index.html").write_text(html)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", type=Path, required=True, help="the benchmark data repository")
    main(parser.parse_args().data / "results_time_curves" / "data")

"""Shared machinery of the benchmark drivers: mesh and fiber conversion to tetgen, the input deck,
the run, probing the displacement, the cohort, RED and the records.

Lengths are in m, as in the benchmark data, so the decks read the meshes with unit 1.
"""
import csv
import pickle
import subprocess
import time
from pathlib import Path
from string import Template

import h5py
import matplotlib
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = Path(__file__).resolve().parent
RESULTS = REPO / "results"
MATERIAL = 30
ENDO, EPI, BASE = 1, 2, 3                  # surface tags of the monoventricle
# The order of the cohort in the paper's figures and tables, with the pickle name of each code.
# The data also holds a `simvascular` pickle, a copy of `simvascular_p2` that the paper does not use.
COHORT = {"CARPentry": "carpentry", "Ambit": "ambit", "4C": "4c", "Simula": "simula", "CHimeRA": "chimera",
          "CHeart": "cheart", "lifex": "lifex", "SimVascular P1": "simvascular_p1p1",
          "SimVascular P2": "simvascular_p2", "COMSOL": "comsol"}
# The cohort of the biventricle, Table 12. Its `simvascular` pickle is SimVascular's P1 result.
BIV_COHORT = {"CARPentry": "carpentry", "Ambit": "ambit", "4C": "4c", "Simula": "simula", "CHimeRA": "chimera",
              "CHeart": "cheart", "lifex": "lifex", "SimVascular P1": "simvascular", "COMSOL": "comsol"}
SAMPLES = 101                              # every 10 ms over [0, 1] s, eq. (21)
COMPONENTS = ("ux", "uy", "uz")
# T10 nodes: the vertices, then the midpoints of these vertex pairs (CBQuadraticTetBasis).
T10_EDGES = ((0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3))
# Barycentric coordinates of the points a T10 bases line holds frames at. The decks need quadrature
# degree 5, for which the T10 element takes a single frame, at the centroid.
BASIS_POINTS = np.array([[0.25] * 4])


def write_tetgen_t4(h5, directory, stem, tags=(ENDO, EPI, BASE)):
    """Write the tetrahedra of a dolfin mesh file and its faces of tags as tetgen stem.node, .ele and
    .sur, with every face normal pointing out of the volume. The faces keep the tag as material and
    surface index; other faces are dropped."""
    with h5py.File(h5, "r") as f:
        points = f["Mesh/mesh/geometry"][:]
        tets = f["Mesh/mesh/topology"][:].astype(int)
        faces = f["MeshFunction/0/mesh/topology"][:].astype(int)
        values = f["MeshFunction/0/values"][:].ravel()
    tagged = np.isin(values, tags)
    faces, values = faces[tagged], values[tagged].astype(int)

    # The vertex of the tetrahedron opposite each tagged face, which lies inside the volume.
    opposite = {}
    for tet in tets:
        for k in range(4):
            opposite[tuple(sorted(np.delete(tet, k)))] = tet[k]
    inside = points[[opposite[tuple(sorted(face))] for face in faces]]
    a, b, c = (points[faces[:, k]] for k in range(3))
    flip = np.einsum("ij,ij->i", np.cross(b - a, c - a), inside - a) > 0
    faces[flip] = faces[flip][:, [0, 2, 1]]

    np.savetxt(directory / f"{stem}.node", np.column_stack([np.arange(1, len(points) + 1), points, np.zeros(len(points))]),
               fmt=["%d", "%.17g", "%.17g", "%.17g", "%d"], header=f"{len(points)} 3 1 0", comments="")
    np.savetxt(directory / f"{stem}.ele", np.column_stack([np.arange(1, len(tets) + 1), tets + 1, np.full(len(tets), MATERIAL)]),
               fmt="%d", header=f"{len(tets)} 4 1", comments="")
    np.savetxt(directory / f"{stem}.sur", np.column_stack([np.arange(1, len(faces) + 1), faces + 1, values, values]),
               fmt="%d", header=f"{len(faces)} 3 2", comments="")


def convert_t10(binary, t4, t10):
    """The T10 mesh of the T4 tetgen files with prefix t4, written with prefix t10 by CardioMechanics'
    ConvertT4toT10. Its mid-edge nodes sit at the edge midpoints."""
    subprocess.run([str(binary), *(f"{t4}.{ext}" for ext in ("node", "ele", "sur")), str(t10)],
                   check=True, stdout=subprocess.DEVNULL)


def write_bases(points, f, s, directory, stem):
    """Write stem.bases for the T10 mesh stem in directory from the degree-2 fiber fields f and s
    given at points, all of shape (n, 3).

    Each node takes the fields at the nearest point, which must be its own. The tolerance admits
    points rounded to six digits, as VTK files of the data store them, far below the node spacing.
    The fields are interpolated with the T10 shape functions to BASIS_POINTS; then f is normalized,
    s orthogonalized against f and normalized, and n = f x s, which makes each frame right-handed. The paper's n, eq. (15), is s x f, the opposite sign, which
    the Holzapfel law does not see."""
    nodes = np.loadtxt(directory / f"{stem}.node", skiprows=1, usecols=(1, 2, 3))
    elements = np.loadtxt(directory / f"{stem}.ele", skiprows=1, usecols=range(1, 11), dtype=int) - 1
    X = nodes[elements]
    for k, (i, j) in enumerate(T10_EDGES, 4):
        assert np.allclose(X[:, k], (X[:, i] + X[:, j]) / 2, rtol=0, atol=1e-12), f"node {k} is not the midpoint of {i}, {j}"

    cloud = vtk.vtkPolyData()
    cloud.SetPoints(vtk.vtkPoints())
    cloud.GetPoints().SetData(numpy_to_vtk(points, deep=True))
    locator = vtk.vtkKdTreePointLocator()
    locator.SetDataSet(cloud)
    locator.BuildLocator()
    index = np.array([locator.FindClosestPoint(x) for x in nodes])
    gap = np.linalg.norm(points[index] - nodes, axis=1).max()
    assert gap < 1e-5, f"a node is {gap} m off the nearest fiber point"
    assert len(np.unique(index)) == len(nodes), "two nodes share a fiber point"

    bary = BASIS_POINTS
    N = np.hstack([bary * (2 * bary - 1), np.stack([4 * bary[:, i] * bary[:, j] for i, j in T10_EDGES], axis=1)])
    f, s = (np.einsum("qa,eaj->eqj", N, v[index][elements]) for v in (f, s))
    f /= np.linalg.norm(f, axis=2, keepdims=True)
    s -= np.sum(s * f, axis=2, keepdims=True) * f
    s /= np.linalg.norm(s, axis=2, keepdims=True)
    frames = np.stack([f, s, np.cross(f, s)], axis=2).reshape(len(elements), -1)
    np.savetxt(directory / f"{stem}.bases", np.column_stack([np.arange(1, len(frames) + 1), frames]),
               fmt=["%d"] + ["%.12g"] * frames.shape[1], header=f"{len(frames)} {len(BASIS_POINTS)}", comments="")


def write_deck(template, wd, **values):
    (wd / "settings.xml").write_text(Template((REPO / "templates" / template).read_text()).substitute(values))


def run(binary, wd, ranks):
    """Run CardioMechanics on wd/settings.xml. Returns the wall time in s and whether the run finished."""
    start = time.perf_counter()
    with open(wd / "stdout.txt", "w") as out:
        proc = subprocess.run(["mpirun", "-np", str(ranks), str(binary), "-settings", "settings.xml"], cwd=wd,
                              stdout=out, stderr=subprocess.STDOUT)
    wall = time.perf_counter() - start
    # A failed simulation can still exit 0, so success is the binary's own report of it.
    finished = "Simulation finished successfully" in (wd / "stdout.txt").read_text(errors="replace")
    return wall, proc.returncode == 0 and finished


def newton_iterations(wd):
    return int(np.loadtxt(wd / "out" / "run_data" / "SnesIterations.dat", ndmin=2)[:, 1].sum())


def last_time(wd):
    """The time of the last step the run wrote, which for a failed run is the step that failed."""
    return float(np.loadtxt(wd / "out" / "run_data" / "SnesIterations.dat", ndmin=2)[-1, 0])


def probe(vtus, node_file, points):
    """Displacement (len(vtus), len(points), 3) at points of the reference configuration.

    CardioMechanics exports the deformed coordinates and the tetgen node of each point as PointID.
    The cell holding each probe point and its shape function weights are found once, in the reference
    configuration rebuilt from the tetgen nodes, and applied to every exported state."""
    def read(vtu):
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(str(vtu))
        reader.Update()
        grid = reader.GetOutput()
        return grid, vtk_to_numpy(grid.GetPoints().GetData()).astype(float)

    grid, _ = read(vtus[0])
    pid = vtk_to_numpy(grid.GetPointData().GetArray("PointID")).astype(int)
    reference = np.loadtxt(node_file, skiprows=1, usecols=(1, 2, 3))[pid]
    grid.GetPoints().SetData(numpy_to_vtk(reference, deep=True))
    locator = vtk.vtkStaticCellLocator()
    locator.SetDataSet(grid)
    locator.BuildLocator()
    stencils = []
    for point in points:
        cell_id = locator.FindCell(point)
        assert cell_id >= 0, f"probe point {point} lies outside the mesh"
        cell = grid.GetCell(cell_id)
        weights = [0.0] * cell.GetNumberOfPoints()
        cell.EvaluatePosition(point, [0.0] * 3, vtk.reference(0), [0.0] * 3, vtk.reference(0.0), weights)
        stencils.append(([cell.GetPointId(i) for i in range(cell.GetNumberOfPoints())], np.asarray(weights)))

    u = np.empty((len(vtus), len(points), 3))
    for t, vtu in enumerate(vtus):
        displacement = read(vtu)[1] - reference
        u[t] = [w @ displacement[ids] for ids, w in stencils]
    return u


def cohort(data, prefix, probes, codes=COHORT):
    """Displacement (SAMPLES, len(probes), 3) of every code of codes, from the pickles in data whose
    name starts with prefix. Samples are matched by index, as the paper does: two codes put their
    first sample at 1 ms instead of 0."""
    files = {p.stem.split("_group_")[1].lower(): p for p in data.glob(f"{prefix}_group_*")}
    curves = {}
    for code, name in codes.items():
        with open(files[name], "rb") as f:
            d = pickle.load(f)["displacement"]
        curves[code] = np.stack([[np.asarray(d[p][c], dtype=float) for c in COMPONENTS] for p in probes]).transpose(2, 0, 1)
        assert curves[code].shape == (SAMPLES, len(probes), 3), f"{code}: {curves[code].shape}"
    return curves


def red(u, mean):
    """RED of eq. (21) per probe point, for u and mean of shape (SAMPLES, probes, 3)."""
    return np.mean(np.linalg.norm(u - mean, axis=2) / np.linalg.norm(mean, axis=2), axis=0)


def plot(path, title, probes, curves, u):
    """One panel per probe point and component: the cohort, thin, and CardioMechanics, bold."""
    t = np.linspace(0, 1, SAMPLES)
    fig, axs = plt.subplots(3, len(probes), sharex=True, figsize=(6 * len(probes), 10), squeeze=False)
    for code, v in curves.items():
        for j in range(len(probes)):
            for i in range(3):
                axs[i, j].plot(t, v[:, j, i], lw=0.8, label=code)
    for j, p in enumerate(probes):
        axs[0, j].set_title(p)
        axs[2, j].set_xlabel("time [s]")
        for i, c in enumerate(COMPONENTS):
            axs[i, j].plot(t, u[:, j, i], color="black", lw=2, label="CardioMechanics")
    for i, c in enumerate(COMPONENTS):
        axs[i, 0].set_ylabel(f"{c} [m]")
    axs[0, -1].legend(fontsize=7)
    fig.suptitle(title)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=120)
    plt.close(fig)


def evaluate(wd, finished, node_file, probes, curves, name, title):
    """The cohort's worst RED and the outcome of the run in wd, as CSV columns.

    A finished run is probed at probes (name: point), scored against the mean of the cohort curves,
    and passes if its RED at every probe point is no worse than the worst cohort code's. Its overlay
    plot and probe trace go to RESULTS as name.png and name_trace.csv. A failed run is recorded with
    the time it failed at and fails."""
    mean = np.mean(list(curves.values()), axis=0)
    threshold = np.max([red(v, mean) for v in curves.values()], axis=0)
    columns = {f"cohort_max_red_{p}": f"{r:.4f}" for p, r in zip(probes, threshold)}
    if not finished:
        columns.update({"pass": False, "failed_at": f"{last_time(wd):g}"})
        print(f"failed at t = {columns['failed_at']} s, see {wd / 'stdout.txt'}")
        return columns
    vtus = sorted((wd / "out" / "run_vtu").glob("run.*.vtu"), key=lambda p: int(p.stem.split(".")[1]))
    assert len(vtus) == SAMPLES, f"{len(vtus)} exported states, expected {SAMPLES}"
    u = probe(vtus, node_file, list(probes.values()))
    score = red(u, mean)
    plot(RESULTS / f"{name}.png", title, list(probes), curves, u)
    np.savetxt(RESULTS / f"{name}_trace.csv", np.column_stack([np.linspace(0, 1, SAMPLES), u.reshape(len(u), -1)]),
               delimiter=",", fmt="%.10g", comments="", header=",".join(["t"] + [f"{p}_{c}" for p in probes for c in COMPONENTS]))
    print("RED", dict(zip(probes, score.round(4))), "cohort max", dict(zip(probes, threshold.round(4))))
    columns.update({f"red_{p}": f"{r:.4f}" for p, r in zip(probes, score)})
    columns["pass"] = bool(np.all(score <= threshold))
    return columns


def dofs(element_file):
    """Unknowns of a T10 mesh: three per node."""
    return 3 * len(np.unique(np.loadtxt(element_file, skiprows=1, usecols=range(1, 11), dtype=int)))


def provenance(binary_dir):
    """Branch and commit of the git checkout holding the binaries."""
    def git(*command):
        return subprocess.run(["git", "-C", str(binary_dir), *command], capture_output=True, text=True, check=True).stdout.strip()
    return {"branch": git("rev-parse", "--abbrev-ref", "HEAD"), "commit": git("rev-parse", "--short", "HEAD")}


def record(path, row):
    """Append row to the CSV at path, rewriting it when the row brings new columns."""
    rows = []
    if path.exists():
        with open(path) as f:
            rows = list(csv.DictReader(f))
    rows.append(row)
    fields = list(dict.fromkeys(k for r in rows for k in r))
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, restval="")
        writer.writeheader()
        writer.writerows(rows)

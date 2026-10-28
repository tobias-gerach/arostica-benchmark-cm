"""Fibers of the monoventricle at the degree-2 nodes, by the benchmark's own toolkit.

Runs inside the toolkit's dolfin Docker image: `python3 dolfin_fibers.py MESH.xdmf OUT.npz`. The
toolkit solves the transmural Laplace problem and interpolates the analytic f, s, n of the paper's
eq. (15) into P2, which evaluates them at the vertices and edge midpoints. Writes those points as X
and the fields f and s, each of shape (points, 3). The bases take n from f and s.
"""
import sys

import numpy as np
from cardiac_benchmark_toolkit.ellipsoid_fiber_generation import build_orientations, transmural_distance_problem
from cardiac_benchmark_toolkit.utils import read_mesh

DEGREE = 2

mesh, boundaries = read_mesh(sys.argv[1])
fields = build_orientations(mesh, transmural_distance_problem(mesh, boundaries, DEGREE), DEGREE)
components = [field.sub(i, deepcopy=True) for field in (fields.fiber, fields.sheet) for i in range(3)]
X = components[0].function_space().tabulate_dof_coordinates()
# Every collapsed component space must list its degrees of freedom at the same points.
assert all(np.array_equal(c.function_space().tabulate_dof_coordinates(), X) for c in components)
values = np.stack([c.vector().get_local() for c in components], axis=1)
np.savez(sys.argv[2], X=X, f=values[:, 0:3], s=values[:, 3:6])

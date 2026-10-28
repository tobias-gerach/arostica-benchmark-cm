# Biventricle: pass/fail

T10, RhoInf 0.5, quadrature degree 5. A run passes if its RED at every probe point is no worse than the worst cohort member's (cohort max), computed from the pickles. The verdict names the probe points that fail. The worst RED of Table 12 is listed next to it.

| Mesh | DOFs | RED p0 | RED p1 | RED p2 | cohort max p0 | cohort max p1 | cohort max p2 | Table 12 p0 | Table 12 p1 | Table 12 p2 | Verdict | Commit | Plot |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| coarse | 227364 | 0.1307 | 0.1249 | 0.3309 | 0.9146 | 0.8484 | 1.4719 | 0.915 | 0.848 | 1.472 | pass | b5cb05e | [png](biv_coarse.png) |
| fine | 574728 | 0.1231 | 0.0904 | 0.3637 | 0.9974 | 0.5624 | 1.8912 | 1.019 | 0.513 | 1.769 | pass | b5cb05e | [png](biv_fine.png) |

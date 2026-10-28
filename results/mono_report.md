# Monoventricle: pass/fail

T10, RhoInf 0.5, quadrature degree 5. A run passes if its RED at every probe point is no worse than the worst cohort member's (cohort max). The verdict names the probe points that fail. Cases A to C list the worst RED of Table 11 next to it.

| Case | DOFs | RED p0 | RED p1 | cohort max p0 | cohort max p1 | Table 11 p0 | Table 11 p1 | Verdict | Commit | Plot |
|---|---|---|---|---|---|---|---|---|---|---|
| 0A | 89679 | 0.1892 | 0.2165 | 1.4742 | 1.5450 | – | – | pass | b5cb05e | [png](mono_0A.png) |
| 0B | 89679 | 0.1488 | 0.1482 | 0.2875 | 0.3702 | – | – | pass | b5cb05e | [png](mono_0B.png) |
| 1 | 89679 | 0.0431 | 0.0669 | 0.1705 | 0.2262 | – | – | pass | b5cb05e | [png](mono_1.png) |
| A | 89679 | 0.0683 | 0.1213 | 0.2757 | 0.3713 | 0.276 | 0.371 | pass | b5cb05e | [png](mono_A.png) |
| B | 89679 | 0.1391 | 0.0987 | 0.3603 | 0.3519 | 0.360 | 0.352 | pass | b5cb05e | [png](mono_B.png) |
| C | 89679 | 0.0426 | 0.0644 | 0.1556 | 0.2111 | 0.156 | 0.211 | pass | bf32b16 | [png](mono_C.png) |

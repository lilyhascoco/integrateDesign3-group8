# Quadratic Convergence Experiments

## Step A: Objective Definition and the Star-Model Decision

The historical clique implementation rebuilds its quadratic system in every
outer iteration. For each pin pair `(i, j)` in a net, it computes separate
directional weights from the current absolute pin locations:

```text
w_ij^x = net_weight / max(|p_i^x - p_j^x|, min_distance)
w_ij^y = net_weight / max(|p_i^y - p_j^y|, min_distance)
```

Thus, for one direction, its reweighted clique objective is:

```text
Cost_clique = sum_nets sum_(i<j in net) w_ij(d_ij) * (q_i - q_j)^2
```

where `q` is the pin coordinate in that direction. The outer loop rebuilds
these weights after each solve. This is an iterative reweighted quadratic
surrogate intended to pull long pin separations together and thereby
approximate wirelength; it is not a fixed-weight quadratic objective.

The current degree >= 3 star implementation instead uses one directional
weight per net, derived from the current net span:

```text
w_net^x = degree * net_weight / max(max(p^x) - min(p^x), min_distance)
w_net^y = degree * net_weight / max(max(p^y) - min(p^y), min_distance)

Cost_star = sum_nets sum_(i in net) w_net * (q_i - q_c)^2
q_c = weighted_average(q_i)
```

Eliminating `q_c` produces a uniform-weight clique for that net. It cannot
reproduce the historical clique when pair weights differ within the same net,
because the historical `w_ij` depends on each individual pair distance while
the star has only one `w_net`.

Decision: current degree >= 3 star placement is an **approximation**, not an
exact transformation of the historical distance-reweighted clique objective.
The prior two-pin sample is not evidence for the star path because two-pin
nets continue to use the original direct edge. Subsequent validation must use
a degree >= 3 net and quantify coordinate, HPWL, and residual differences.

## Stage 0: Clique / Star Equivalence Gate

Date: 2026-10-24

Input: `benchmarks/sample/sample.aux`

Commands:

```powershell
# Historical clique baseline at HEAD (commit 56b51cc)
placer.exe benchmarks/sample/sample.aux --method quadratic --no-sa --output e:\clique-sample.pl

# Current star-model worktree
placer.exe benchmarks/sample/sample.aux --method quadratic --no-sa --output e:\star-sample.pl
```

Parameters:

| Version | Net model | Outer iterations | BiCGSTAB limit | Residual definition |
| --- | --- | ---: | ---: | --- |
| Historical baseline | all pin-pair clique | 5 used | 500 | preconditioned residual infinity norm |
| Current worktree | degree >= 3 star, degree 2 direct edge | 5 used | 50 | preconditioned residual infinity norm |

Results:

| Metric | Clique | Current model |
| --- | ---: | ---: |
| Final HPWL | 20.19 | 20.19 |
| Final X residual | 8.72e-06 | 8.72e-06 |
| Final Y residual | 9.01e-13 | 9.01e-13 |
| Runtime | 0.033 ms | 0.068 ms |
| Bin-overlap overflow | 95.31% | 95.31% |
| Maximum absolute cell-coordinate difference | 0 | 0 |

The exported `.pl` files had identical SHA-256 hashes:

```text
42EDF510686060A4FA91272AF476B9FBCDC8AA86003D03CCFAF812C65D4BFA83
```

Conclusion: the sample contains only two-pin nets, whose current direct-edge
path is exactly identical to the historical clique formulation. The Stage 0
sample gate passes. This does not by itself prove equivalence for degree >= 3
nets; that limitation must be considered during Stage 1 diagnosis.

# Quadratic Convergence Experiments

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

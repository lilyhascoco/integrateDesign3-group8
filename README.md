# Placement Course Project

一个面向芯片物理设计 Placement 课程的 C++ 工程骨架，当前阶段重点实现：

- BookShelf 输入解析
- 统一内部数据模型 `PlacementDB + PlacementState`
- 基础统计与初始 HPWL 计算

## Build

```bash
cmake -S . -B build
cmake --build build
```

## Run

```bash
./build/placer <path-to-benchmark.aux>
./build/placer <path-to-benchmark.aux> --check-only
```

## Task 4: Initial placement

See [算法与运行说明](docs/任务四-聚类驱动布局.md).

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
ctest --test-dir build --output-on-failure
# If Python was not discovered by CMake:
python3 tests/cli_integration.py build/placer
./build/placer benchmarks/sample/sample.aux --method compare --output build/sample.pl
./build/placer benchmarks/adaptec1/adaptec1.aux --method compare --no-sa
```

`--method random|cluster|quadratic|compare` selects an initializer or runs all three
independently from the same input. Default: random, without simulated annealing.
Use `--sa` to run SA separately after each selected initializer; `--no-sa` disables
only SA. `--check-only` parses the input without placing cells.
Compare exports are named `sample.random.pl`, `sample.cluster.pl`, and
`sample.quadratic.pl`. Exports contain initial placement, before optional SA.
All methods use the same center-offset HPWL evaluator. Bin overflow is a density
proxy, not exact overlap or a legality guarantee. Quadratic residual warnings must
be reviewed independently of boundary checks.

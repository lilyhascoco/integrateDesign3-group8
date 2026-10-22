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

## Task 4: Clustering-driven placement

See [算法与 VS Code 运行说明](docs/任务四-聚类驱动布局.md) and [实测结果](docs/实验结果.md).

```bash
./build/placer benchmarks/adaptec1/adaptec1.aux --method compare --no-sa
```

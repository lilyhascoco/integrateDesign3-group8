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

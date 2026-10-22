# Placement Course Project

一个面向芯片物理设计 Placement 课程的 C++ 工程，覆盖：

- BookShelf 输入解析（`.aux/.nodes/.nets/.pl/.scl/.wts`）
- 统一内部数据模型 `PlacementDB + PlacementState`
- 基础统计与 HPWL 计算
- 全局布局（任务五，ePlace 风格静电学模型）
  - Filler 初始化与 m×m 密度 Bin 网格
  - 密度 DCT 变换、电势/电场求解与密度梯度
  - WA 线长平滑模型与线长梯度
  - Nesterov 加速梯度法（Lipschitz 步长预测）
  - λ 自适应更新与密度溢出率停止准则
- 局部模拟退火（可选的后续优化阶段）

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

常用选项：

```bash
--iterations <N>        最大全局布局迭代次数（默认 200）
--target-density <v>    目标密度（默认 0.9）
--target-overflow <v>   密度溢出率停止阈值（默认 0.1）
--bins <N>              Bin 网格边长上限（默认 512）
--seed <N>              随机种子（默认 2024）
--output <file.pl>      导出全局布局结果（UCLA .pl 格式）
--local-sa              全局布局后继续运行局部模拟退火
--quiet                 关闭逐轮日志
```

例如：

```bash
./build/placer benchmarks/sample/sample.aux --output build/sample.gp.pl
./build/placer benchmarks/adaptec1/adaptec1.aux --iterations 200 --quiet
```

## 任务五：全局布局（静电学模型）

目标函数与梯度：

```
f(x, y) = HPWL(x, y) + lambda * D(x, y)
grad f  = grad HPWL + lambda * grad D
```

实现流程（`src/algorithms/global_placement.cpp`）：

1. **FillerInit**：filler 总面积 = targetDensity·coreArea − (可移动面积 + 终端在核心区面积)，
   单个 filler 边长取中间 80% 单元的平均面积开方。
2. **BinInit**：m×m 网格，m 取不小于 sqrt(节点数) 的 2 的幂，上限 512。
3. **密度图**：标准单元贡献 overlapArea；Macro/Filler/Terminal 再乘 targetDensity；
   小于 bin 的节点按 bin 尺寸窗口均匀涂抹（质量守恒），行间空白作为 dark area 计入。
4. **电势与电场**：密度图经 DCT-II 变换得到系数 a(j,k)，用
   `a(j,k) / (w_j^2 + w_k^2)`（w = pi·k/m）做 DCT-III 得到电势 psi；
   电场取中心差分 −∇psi，节点密度梯度按节点覆盖 bin 的重叠面积加权平均。
5. **WA 线长模型**：`gamma(tau) = 8·w_b·10^((20/9)(tau−0.1)−1)`，用 b±/c±/t± 计算平滑线长梯度。
6. **λ 更新**：首轮 `lambda0 = Σ|grad HPWL| / Σ|grad D|`，之后
   `lambda_k = mu·lambda_{k−1}`，`mu = 1.1^(1 − ΔHPWL/3.5e5)`，截断到 [0.75, 1.1]。
7. **Nesterov 加速梯度法**：u/v 双序列 + Lipschitz 常数预测步长
   `L = ||grad f(v_k) − grad f(v_{k−1})|| / ||v_k − v_{k−1}||`，并对单步位移做上限保护。
8. **停止准则**：密度溢出率 `tau = 总溢出面积 / 总单元面积` < targetOverflow（默认 0.1）。

## 验证结果

| benchmark | Bin 网格 | filler | 迭代 | 初始 HPWL | 最终 HPWL | 初始溢出 | 最终溢出 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| adaptec1 | 512×512 | 96425 | 137 | 2.2130e9 | 2.0573e9 | 45.05% | 9.93% |

sample 是一个只有 2 个可移动单元、核心区高度等于单元高度的退化算例（只有 x 方向可动，
且固定终端占据两个 bin 列），其离散化后可达的最小溢出率本身高于默认 10%，
因此该算例只用于功能冒烟测试：

```bash
./build/placer benchmarks/sample/sample.aux --iterations 200 --target-overflow 0.3
```


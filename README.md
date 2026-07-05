# Engineering Modeling and Simulation MATLAB Cases

这是我在“工程问题建模与仿真”课程中整理出来的两个 MATLAB case。这个仓库不是完整作业归档，而是一个公开展示版本：保留能说明建模思路、仿真方法、运行入口和结果的代码与图表，去掉课程题面、提交报告、参考资料和过程性备份。

## Cases

| Case | Focus | Main MATLAB Work |
| --- | --- | --- |
| `case1` | 测量点方案优化 | 三次样条插值、代价函数评估、遗传算法/局部搜索、结果可视化 |
| `case2` | 系统可靠性仿真 | Monte Carlo 可靠度/MTTF 估计、瞬时可用度枚举、失效后 revival 诊断 |

## Repository Layout

```text
case1/
  code/                       MATLAB source and input CSV data
  results/submit_example/     selected output tables and figures

case2/
  code/main_case2_submit.m    standalone MATLAB script
  results/                    selected numerical results
  figures/                    selected output plots

docs/PUBLICATION_NOTES.md     what was included/excluded for public release
```

## How to Run

The source environment was MATLAB R2024b. The scripts use standard MATLAB language features and built-in plotting/statistical utilities; no third-party toolbox is required by the public entry points.

### Case 1

```matlab
cd case1/code
run_case1_submit
```

`run_case1_submit` uses the public verification profile from the original submission package. It regenerates a timestamped output folder under `case1/code/results/`. This is a real search workflow and can take a while; `case1/results/submit_example/` is kept so the output can be inspected without rerunning the full search.

### Case 2

```matlab
cd case2/code
main_case2_submit
```

The default run uses `S = 100000` Monte Carlo samples. For a quick smoke run:

```matlab
cd case2/code
setenv("CASE2_SAMPLE_SIZE", "2000")
setenv("CASE2_OUTPUT_TAG", "smoke")
main_case2_submit
```

## Selected Results

Case 1 public example:

- Best point count in the quick verification run: `7`
- Selected temperatures: `[-17, -4, 8, 26, 36, 54, 68]`
- `C_train = 524.112`, `C_testA = 524.875`
- Train/TestA mean absolute error: `0.2269 / 0.2279`

Case 2 public example:

- Reliability target: `w = 25000 hours`
- Monte Carlo sample size: `100000`
- Best reliability node count: `n = 9`, `R(25000) = 0.785350`
- Best MTTF point estimate: `n = 10`, `MTTF = 60304.81 hours`
- The public conclusion keeps `n = 9` as the reliability-first recommendation.

## Public Boundary

This repository intentionally excludes course PDFs, problem statements, full reports, local absolute paths, submission archives, historical backup folders and reference materials. The goal is to show the MATLAB modeling/simulation work itself, not to mirror the whole course folder.

## Version

`v1.0-public` collects the two verified coursework cases into a compact GitHub-ready form.

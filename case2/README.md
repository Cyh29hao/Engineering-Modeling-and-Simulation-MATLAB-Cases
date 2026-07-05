# Case 2: Reliability and MTTF Simulation

This case models a multi-node system reliability problem. The public MATLAB script estimates reliability and MTTF with Monte Carlo simulation, compares reliability against an instantaneous-state availability approximation, and records whether a path can become workable again after a first failure.

## Files

- `code/main_case2_submit.m` is a standalone MATLAB entry point.
- `results/case2_submit_results.csv` contains the selected numerical result table.
- `results/case2_submit_summary.txt` records the key run parameters and conclusions.
- `figures/` contains the reliability and MTTF plots generated from the selected run.

## Run

```matlab
cd case2/code
main_case2_submit
```

The default sample size is `100000`. For a quicker check:

```matlab
cd case2/code
setenv("CASE2_SAMPLE_SIZE", "2000")
setenv("CASE2_OUTPUT_TAG", "smoke")
main_case2_submit
```

## Public Example Result

- reliability target: `w = 25000 hours`
- default sample size: `100000`
- best reliability node count: `n = 9`
- `R(25000) = 0.785350`
- best MTTF point estimate: `n = 10`
- `MTTF = 60304.81 hours`

The public recommendation keeps `n = 9` because the task prioritizes reliability at the target time, while the MTTF advantage at `n = 10` is relatively small in the recorded run.

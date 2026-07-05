# Case 1: Measurement Scheme Optimization

This case searches for a compact set of measurement temperatures. The public code evaluates candidate schemes with cubic spline interpolation and optimizes the point set with a genetic-search style workflow plus local refinement.

## Files

- `code/run_case1_submit.m` is the public entry point.
- `code/case1_submit_config.m` defines the search profiles, including `submit` and `smoke`.
- `code/case1_run_search.m` coordinates loading data, running search, writing tables and exporting figures.
- `code/dataform_train2026.csv` and `code/dataform_testA2026.csv` are the input data used by the public script.
- `results/submit_example/` keeps one selected output run for quick inspection.

## Run

```matlab
cd case1/code
run_case1_submit
```

The script creates a new timestamped directory under `case1/code/results/`. It is a real search workflow, so runtime depends on the local MATLAB environment and selected profile. The public repository keeps one selected output run under `results/submit_example/` for quick inspection.

## Public Example Result

The included example output uses the public verification profile:

- point-count range: `4` to `12`
- population per island: `20`
- max generations: `5`
- best point count: `7`
- selected temperatures: `[-17, -4, 8, 26, 36, 54, 68]`
- `C_train = 524.112`
- `C_testA = 524.875`

This is presented as a reproducible coursework verification sample, not as a claim of global optimality.

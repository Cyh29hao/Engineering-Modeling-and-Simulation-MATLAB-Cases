function run_case1_submit()
%RUN_CASE1_SUBMIT Quick verification entry point for Case1 submission.
% The submitted verification code uses population size 20 and 5 generations.

scriptPath = mfilename('fullpath');
if isempty(scriptPath)
    caseDir = pwd;
else
    caseDir = fileparts(scriptPath);
end

% MATLAB may keep a stale function cache immediately after a zip is extracted
% or copied into a Chinese-path folder. Refreshing here avoids "function not
% recognized" errors when the teacher runs the demo from the GUI.
cd(caseDir);
addpath(caseDir, '-begin');
rehash;

requiredFiles = {'case1_submit_config.m', 'case1_run_search.m', ...
    'case1_load_data.m', 'case1_evaluate_scheme.m', 'case1_evaluate_spline.m', ...
    'case1_ga_optimize.m'};
missingFiles = requiredFiles(~cellfun(@(name) isfile(fullfile(caseDir, name)), requiredFiles));
if ~isempty(missingFiles)
    error('Case1:MissingFiles', 'Missing required file(s): %s', strjoin(missingFiles, ', '));
end

if exist('case1_submit_config', 'file') ~= 2
    rehash path;
end
if exist('case1_submit_config', 'file') ~= 2
    error(['MATLAB cannot find case1_submit_config.m. Please make the extracted ', ...
        'package folder the current folder, then run: rehash; clear functions; run_case1_submit']);
end

cfg = case1_submit_config('submit');
totalPopulation = cfg.search.populationSizePerIsland * cfg.search.islandCount;

fprintf('===== Case1 submission demo: required quick verification =====\n');
fprintf('Required GA settings: population=20, generations=5\n');
fprintf('Actual GA settings: population=%d, generations=%d\n', ...
    totalPopulation, cfg.search.maxGenerations);
fprintf('Interpolation method = MATLAB interp1 with spline/extrap\n\n');
fprintf(['Submit profile note: this is the required 20-by-5 verification demo; ', ...
    'it keeps only light local adjustment and disables the final heavy refinements.\n']);
fprintf(['Formal report profile note: main PDF results use the fast profile ', ...
    '(3 islands x 20 individuals x 18 generations + archive refinement).\n\n']);

summary = case1_run_search(cfg);

data = case1_load_data(caseDir);
baselineScheme = zeros(1, data.train.pointCount);
baselineScheme([4, 19, 38, 54, 79, 88]) = 1;

detailOptions.costQ = cfg.eval.costQ;
detailOptions.needDetail = true;
detailOptions.enablePruning = false;
detailOptions.method = cfg.eval.method;

baselineTrain = case1_evaluate_scheme(data.train, baselineScheme, detailOptions);
baselineTestA = case1_evaluate_scheme(data.testA, baselineScheme, detailOptions);

fprintf('\n===== Baseline spline result =====\n');
fprintf('baseline indices = [4, 19, 38, 54, 79, 88]\n');
fprintf('baseline temperatures = %s\n', mat2str(baselineTrain.selectedTemperature));
fprintf('baseline C_train = %.3f\n', baselineTrain.meanCost);
fprintf('baseline C_testA = %.3f\n', baselineTestA.meanCost);

fprintf('\n===== Quick GA result =====\n');
fprintf('best selected temperatures = %s\n', mat2str(summary.bestTrainEval.selectedTemperature));
fprintf('quick C_train = %.3f\n', summary.bestTrainEval.meanCost);
fprintf('quick C_testA = %.3f\n', summary.bestTestEval.meanCost);
fprintf('quick gap = %.3f\n', abs(summary.bestTrainEval.meanCost - summary.bestTestEval.meanCost));
fprintf(['Interpretation: the quick run proves the code path and trend plot are reproducible; ', ...
    'the PDF uses the longer fast run for the reported optimum.\n']);

fprintf('\n===== Demo output files =====\n');
fprintf('output directory = %s\n', summary.resultDir);
fprintf('GA trend figure = %s\n', fullfile(summary.resultDir, 'ga_trend_min_mean.png'));
fprintf('search history CSV = %s\n', fullfile(summary.resultDir, 'search_history.csv'));
fprintf('summary text = %s\n', fullfile(summary.resultDir, 'case1_submit_summary.txt'));
end

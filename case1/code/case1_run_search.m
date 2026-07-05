function resultSummary = case1_run_search(cfg)
%CASE1_RUN_SEARCH Run the clean Case1 submit verification search.

if nargin < 1 || isempty(cfg)
    cfg = case1_submit_config('balanced');
elseif ischar(cfg) || isstring(cfg)
    cfg = case1_submit_config(char(cfg));
end

clc;
close all;

caseDir = fileparts(mfilename('fullpath'));
addpath(caseDir);

runTag = datestr(now, 'yyyymmdd_HHMMSS');
resultDir = fullfile(caseDir, 'results', [runTag, '_', cfg.meta.profile]);
if ~exist(resultDir, 'dir')
    mkdir(resultDir);
end

data = case1_load_data(caseDir);

fprintf('===== Case1 clean submission search (%s) =====\n', cfg.meta.profile);
fprintf('train samples = %d, testA samples = %d, grid points = %d\n', ...
    data.train.sampleCount, data.testA.sampleCount, data.train.pointCount);
fprintf('search point-count range = [%d, %d]\n\n', ...
    cfg.search.minPointCount, cfg.search.maxPointCount);

optimizerResult = case1_ga_optimize(data.train, cfg, resultDir);

detailOptions.costQ = cfg.eval.costQ;
detailOptions.needDetail = true;
detailOptions.enablePruning = false;
detailOptions.method = cfg.eval.method;

bestTrainEval = case1_evaluate_scheme(data.train, optimizerResult.bestScheme, detailOptions);
bestTestEval = case1_evaluate_scheme(data.testA, optimizerResult.bestScheme, detailOptions);

countResults = evaluate_count_archive(optimizerResult.bestByCount, data, cfg);

write_result_tables_SUBMIT(optimizerResult, bestTrainEval, bestTestEval, countResults, resultDir);
write_summary_text_SUBMIT(cfg, data, optimizerResult, bestTrainEval, bestTestEval, countResults, resultDir);
create_result_figures_SUBMIT(optimizerResult, countResults, resultDir);

resultSummary.resultDir = resultDir;
resultSummary.optimizerResult = optimizerResult;
resultSummary.bestTrainEval = bestTrainEval;
resultSummary.bestTestEval = bestTestEval;
resultSummary.countResults = countResults;

fprintf('===== Case1 submit final best result =====\n');
fprintf('best point count = %d\n', bestTrainEval.selectedCount);
fprintf('selected indices = %s\n', vector_to_string_SUBMIT(bestTrainEval.selectedIdx));
fprintf('selected temperatures = %s\n', vector_to_string_SUBMIT(bestTrainEval.selectedTemperature));
fprintf('C_train = %.3f\n', bestTrainEval.meanCost);
fprintf('C_testA = %.3f\n', bestTestEval.meanCost);
fprintf('gap = %.3f\n', abs(bestTrainEval.meanCost - bestTestEval.meanCost));
fprintf('Outputs written to:\n%s\n', resultDir);
end

function countResults = evaluate_count_archive(bestByCount, data, cfg)
detailOptions.costQ = cfg.eval.costQ;
detailOptions.needDetail = true;
detailOptions.enablePruning = false;
detailOptions.method = cfg.eval.method;

validMask = arrayfun(@(x) ~isempty(x.scheme), bestByCount);
validRecords = bestByCount(validMask);
countResults = repmat(struct(), numel(validRecords), 1);

for i = 1:numel(validRecords)
    trainEval = case1_evaluate_scheme(data.train, validRecords(i).scheme, detailOptions);
    testEval = case1_evaluate_scheme(data.testA, validRecords(i).scheme, detailOptions);

    countResults(i).pointCount = trainEval.selectedCount;
    countResults(i).scheme = validRecords(i).scheme;
    countResults(i).trainEval = trainEval;
    countResults(i).testEval = testEval;
    countResults(i).gap = abs(trainEval.meanCost - testEval.meanCost);
end

if ~isempty(countResults)
    [~, order] = sort([countResults.pointCount], 'ascend');
    countResults = countResults(order);
end
end

function write_result_tables_SUBMIT(optimizerResult, bestTrainEval, bestTestEval, countResults, resultDir)
countCell = cell(numel(countResults) + 1, 10);
countCell(1, :) = {'point_count', 'train_cost', 'testA_cost', 'gap', ...
    'train_measurement_cost', 'train_error_cost', 'train_mean_abs_error', ...
    'testA_mean_abs_error', 'train_max_abs_error', 'selected_temperatures'};

for i = 1:numel(countResults)
    countCell(i + 1, :) = { ...
        countResults(i).pointCount, ...
        countResults(i).trainEval.meanCost, ...
        countResults(i).testEval.meanCost, ...
        countResults(i).gap, ...
        countResults(i).trainEval.meanMeasurementCost, ...
        countResults(i).trainEval.meanErrorCost, ...
        countResults(i).trainEval.meanAbsError, ...
        countResults(i).testEval.meanAbsError, ...
        countResults(i).trainEval.maxAbsError, ...
        vector_to_string_SUBMIT(countResults(i).trainEval.selectedTemperature)};
end
writecell(countCell, fullfile(resultDir, 'count_archive_results.csv'));

pointCell = cell(bestTrainEval.selectedCount + 1, 3);
pointCell(1, :) = {'index', 'temperature', 'binary_gene'};
for i = 1:bestTrainEval.selectedCount
    pointCell(i + 1, :) = { ...
        bestTrainEval.selectedIdx(i), ...
        bestTrainEval.selectedTemperature(i), ...
        optimizerResult.bestScheme(bestTrainEval.selectedIdx(i))};
end
writecell(pointCell, fullfile(resultDir, 'best_scheme_points.csv'));

history = optimizerResult.historyTable;
historyCell = cell(height(history) + 1, width(history));
historyCell(1, :) = history.Properties.VariableNames;
for row = 1:height(history)
    for col = 1:width(history)
        historyCell{row + 1, col} = history{row, col};
    end
end
writecell(historyCell, fullfile(resultDir, 'search_history.csv'));

bestCell = {
    'metric', 'train', 'testA';
    'mean_cost', bestTrainEval.meanCost, bestTestEval.meanCost;
    'mean_measurement_cost', bestTrainEval.meanMeasurementCost, bestTestEval.meanMeasurementCost;
    'mean_error_cost', bestTrainEval.meanErrorCost, bestTestEval.meanErrorCost;
    'mean_abs_error', bestTrainEval.meanAbsError, bestTestEval.meanAbsError;
    'max_abs_error', bestTrainEval.maxAbsError, bestTestEval.maxAbsError;
    'avg_count_over_0p4', bestTrainEval.avgCountOver04, bestTestEval.avgCountOver04;
    'avg_count_over_0p8', bestTrainEval.avgCountOver08, bestTestEval.avgCountOver08;
    'avg_count_over_1p2', bestTrainEval.avgCountOver12, bestTestEval.avgCountOver12;
    'total_count_over_2p0', bestTrainEval.totalCountOver20, bestTestEval.totalCountOver20};
writecell(bestCell, fullfile(resultDir, 'best_scheme_metrics.csv'));
end

function write_summary_text_SUBMIT(cfg, data, optimizerResult, bestTrainEval, bestTestEval, countResults, resultDir)
summaryPath = fullfile(resultDir, 'case1_submit_summary.txt');
fid = fopen(summaryPath, 'w', 'n', 'UTF-8');

fprintf(fid, 'Case1 submit summary\n\n');
fprintf(fid, '1. Profile\n');
fprintf(fid, '- profile: %s\n', cfg.meta.profile);
fprintf(fid, '- point-count range: [%d, %d]\n', cfg.search.minPointCount, cfg.search.maxPointCount);
fprintf(fid, '- restarts: %d\n', cfg.random.restartCount);
fprintf(fid, '- islands: %d\n', cfg.search.islandCount);
fprintf(fid, '- population per island: %d\n', cfg.search.populationSizePerIsland);
fprintf(fid, '- max generations: %d\n\n', cfg.search.maxGenerations);

fprintf(fid, '2. Dataset\n');
fprintf(fid, '- train samples: %d\n', data.train.sampleCount);
fprintf(fid, '- testA samples: %d\n', data.testA.sampleCount);
fprintf(fid, '- grid points: %d\n\n', data.train.pointCount);

fprintf(fid, '3. Best-by-count archive\n');
for i = 1:numel(countResults)
    fprintf(fid, '- k = %d, train = %.3f, testA = %.3f, gap = %.3f, T = %s\n', ...
        countResults(i).pointCount, countResults(i).trainEval.meanCost, ...
        countResults(i).testEval.meanCost, countResults(i).gap, ...
        vector_to_string_SUBMIT(countResults(i).trainEval.selectedTemperature));
end
fprintf(fid, '\n');

fprintf(fid, '4. Final best scheme\n');
fprintf(fid, '- best point count: %d\n', bestTrainEval.selectedCount);
fprintf(fid, '- selected indices: %s\n', vector_to_string_SUBMIT(bestTrainEval.selectedIdx));
fprintf(fid, '- selected temperatures: %s\n', vector_to_string_SUBMIT(bestTrainEval.selectedTemperature));
fprintf(fid, '- C_train: %.3f\n', bestTrainEval.meanCost);
fprintf(fid, '- C_testA: %.3f\n', bestTestEval.meanCost);
fprintf(fid, '- gap: %.3f\n', abs(bestTrainEval.meanCost - bestTestEval.meanCost));
fprintf(fid, '- train measurement cost: %.3f\n', bestTrainEval.meanMeasurementCost);
fprintf(fid, '- train error cost: %.3f\n', bestTrainEval.meanErrorCost);
fprintf(fid, '- train mean abs error: %.4f\n', bestTrainEval.meanAbsError);
fprintf(fid, '- testA mean abs error: %.4f\n', bestTestEval.meanAbsError);
fprintf(fid, '- train max abs error: %.4f\n', bestTrainEval.maxAbsError);
fprintf(fid, '- testA max abs error: %.4f\n', bestTestEval.maxAbsError);
fprintf(fid, '- train total count(|err|>2.0): %.0f\n', bestTrainEval.totalCountOver20);
fprintf(fid, '- testA total count(|err|>2.0): %.0f\n\n', bestTestEval.totalCountOver20);

fprintf(fid, '5. Search outputs\n');
fprintf(fid, '- count_archive_results.csv\n');
fprintf(fid, '- best_scheme_points.csv\n');
fprintf(fid, '- best_scheme_metrics.csv\n');
fprintf(fid, '- search_history.csv\n');
fprintf(fid, '- count_cost_comparison.png\n');
fprintf(fid, '- search_history.png\n');
fprintf(fid, '- ga_trend_min_mean.png\n');

fclose(fid);
end

function create_result_figures_SUBMIT(optimizerResult, countResults, resultDir)
if ~isempty(countResults)
    pointCounts = [countResults.pointCount];
    trainCost = arrayfun(@(x) x.trainEval.meanCost, countResults);
    testCost = arrayfun(@(x) x.testEval.meanCost, countResults);

    fig1 = figure('Visible', 'off');
    plot(pointCounts, trainCost, 'o-r', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
    hold on;
    plot(pointCounts, testCost, 's-b', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
    xlabel('Calibration point count');
    ylabel('Average total cost');
    title('Best archive cost by point count');
    legend('Train', 'TestA', 'Location', 'best');
    grid on;
    saveas(fig1, fullfile(resultDir, 'count_cost_comparison.png'));
    close(fig1);
end

history = optimizerResult.historyTable;
if ~isempty(history)
    fig2 = figure('Visible', 'off');
    yyaxis left;
    plot(history.Iteration, history.GenerationBestSearchCost, 'Color', [0.75 0.75 0.75], 'LineWidth', 1.0);
    hold on;
    plot(history.Iteration, history.GenerationBestExactCost, 'b--', 'LineWidth', 1.2);
    plot(history.Iteration, history.GlobalBestExactCost, 'r-', 'LineWidth', 1.5);
    ylabel('Cost');
    yyaxis right;
    plot(history.Iteration, history.GlobalBestCount, 'k-.', 'LineWidth', 1.1);
    ylabel('Best point count');
    xlabel('Iteration');
    title('Case1 submit search history');
    legend('Generation best search', 'Generation best exact', 'Global best exact', 'Global best count', ...
        'Location', 'best');
    grid on;
    saveas(fig2, fullfile(resultDir, 'search_history.png'));
    close(fig2);

    if ismember('GenerationMeanSearchCost', history.Properties.VariableNames)
        meanSearchCost = history.GenerationMeanSearchCost;
    else
        meanSearchCost = history.GenerationBestSearchCost;
    end

    fig3 = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 980, 680]);
    iter = history.Iteration;
    bestCost = history.GenerationBestSearchCost;
    hasLargeFirstMean = numel(iter) >= 3 && meanSearchCost(1) > 3 * max(meanSearchCost(2:end));

    if hasLargeFirstMean
        axTop = subplot(2, 1, 1, 'Parent', fig3);
        detailIdx = 2:numel(iter);
        plot(axTop, iter(detailIdx), meanSearchCost(detailIdx), 's--', ...
            'Color', [0.85 0.25 0.12], 'LineWidth', 2.0, 'MarkerSize', 7, ...
            'MarkerFaceColor', [0.85 0.25 0.12]);
        hold(axTop, 'on');
        topLimits = expand_limits_SUBMIT(meanSearchCost(detailIdx));
        ylim(axTop, topLimits);
        plot(axTop, [iter(1), iter(1)], topLimits, ':', 'Color', [0.85 0.25 0.12], 'LineWidth', 2.0);
        text(axTop, iter(1) + 0.05, topLimits(2) - 0.12 * diff(topLimits), ...
            sprintf('第1代平均 %.2g，仅作截断标记', meanSearchCost(1)), ...
            'Color', [0.75 0.18 0.08], 'FontWeight', 'bold', 'FontSize', 10);
        ylabel(axTop, '每代平均搜索成本');
        title(axTop, '提交验证 GA 趋势（第1代平均成本截断）');
        legend(axTop, '每代平均成本（第2代起）', '第1代平均成本标记', ...
            'Location', 'southoutside', 'Orientation', 'horizontal');
        grid(axTop, 'on');
        xlim(axTop, [iter(1), iter(end)]);
        set(axTop, 'XTick', iter);

        axBottom = subplot(2, 1, 2, 'Parent', fig3);
        plot(axBottom, iter, bestCost, 'o-', ...
            'Color', [0.10 0.35 0.80], 'LineWidth', 2.0, 'MarkerSize', 7, ...
            'MarkerFaceColor', [0.10 0.35 0.80]);
        ylim(axBottom, expand_limits_SUBMIT(bestCost));
        xlabel(axBottom, '代数');
        ylabel(axBottom, '每代最小搜索成本');
        legend(axBottom, '每代最小成本', 'Location', 'best');
        grid(axBottom, 'on');
        xlim(axBottom, [iter(1), iter(end)]);
        set(axBottom, 'XTick', iter);
    else
        plot(iter, bestCost, 'o-', ...
            'Color', [0.10 0.35 0.80], 'LineWidth', 2.0, 'MarkerSize', 7, ...
            'MarkerFaceColor', [0.10 0.35 0.80]);
        hold on;
        plot(iter, meanSearchCost, 's--', ...
            'Color', [0.85 0.25 0.12], 'LineWidth', 2.0, 'MarkerSize', 7, ...
            'MarkerFaceColor', [0.85 0.25 0.12]);
        xlabel('代数');
        ylabel('搜索成本');
        title('遗传算法成本趋势');
        legend('每代最小成本', '每代平均成本', 'Location', 'best');
        grid on;
        xlim([iter(1), iter(end)]);
        set(gca, 'XTick', iter);
    end
    saveas(fig3, fullfile(resultDir, 'ga_trend_min_mean.png'));
    close(fig3);
end
end

function limits = expand_limits_SUBMIT(values)
values = values(:);
values = values(isfinite(values));
if isempty(values)
    limits = [0, 1];
    return;
end
lo = min(values);
hi = max(values);
if abs(hi - lo) < 1e-9
    pad = max(1, abs(hi) * 0.02);
else
    pad = max(1, 0.12 * (hi - lo));
end
limits = [lo - pad, hi + pad];
end

function out = vector_to_string_SUBMIT(vec)
if isempty(vec)
    out = '[]';
    return;
end

if all(abs(vec - round(vec)) < 1e-9)
    cellText = arrayfun(@(x) sprintf('%d', round(x)), vec, 'UniformOutput', false);
else
    cellText = arrayfun(@(x) sprintf('%.3f', x), vec, 'UniformOutput', false);
end
out = ['[', strjoin(cellText, ', '), ']'];
end


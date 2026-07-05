clear;
clc;
close all;

params.k = 6;
params.nValues = 6:20;
params.w = 25000;
params.maxLife = 250000;
params.meanA = 5.5e4;
params.meanB = 2.2e5;
params.pAType = [0.60, 0.15, 0.25];
params.pBType = [0.65, 0.35];
params.sampleSize = 1e5;
params.seed = 20260519;

sampleOverride = str2double(getenv("CASE2_SAMPLE_SIZE"));
if ~isnan(sampleOverride) && sampleOverride > 0
    params.sampleSize = round(sampleOverride);
end

outputTag = strtrim(getenv("CASE2_OUTPUT_TAG"));
if outputTag == ""
    fileStem = "case2_submit";
else
    fileStem = "case2_submit_" + outputTag;
end

scriptDir = fileparts(mfilename("fullpath"));
if scriptDir == ""
    scriptDir = pwd;
end
submitDir = fileparts(scriptDir);
resultDir = fullfile(submitDir, "results");
figureDir = fullfile(submitDir, "figures");
if ~exist(resultDir, "dir")
    mkdir(resultDir);
end
if ~exist(figureDir, "dir")
    mkdir(figureDir);
end

rng(params.seed, "twister");
runTimer = tic;

nCount = numel(params.nValues);
reliabilityMC = zeros(nCount, 1);
reliabilitySE = zeros(nCount, 1);
mttfMC = zeros(nCount, 1);
mttfSE = zeros(nCount, 1);
availabilityTheory = zeros(nCount, 1);
enumTotalProbability = zeros(nCount, 1);
firstFailureRate = zeros(nCount, 1);
revivalCandidateRate = zeros(nCount, 1);
revivalCandidateGivenFailure = zeros(nCount, 1);
revivalObservedRate = zeros(nCount, 1);
revivalObservedGivenFailure = zeros(nCount, 1);
revivalObservedSE = zeros(nCount, 1);
meanRevivalDelay = nan(nCount, 1);

fprintf("Case 2 sonar reliability simulation\n");
fprintf("S = %d, w = %.0f hours, cap = %.0f hours, seed = %d\n", ...
    params.sampleSize, params.w, params.maxLife, params.seed);
fprintf(" n      R_MC       SE_R        MTTF       SE_MTTF     A_theory   revive|fail enum_err\n");

for idx = 1:nCount
    n = params.nValues(idx);
    [lifetimes, revivalInfo] = simulateCase2Lifetimes(n, params);

    reliabilityMC(idx) = mean(lifetimes >= params.w);
    reliabilitySE(idx) = sqrt(reliabilityMC(idx) * (1 - reliabilityMC(idx)) / params.sampleSize);
    mttfMC(idx) = mean(lifetimes);
    mttfSE(idx) = std(lifetimes, 0) / sqrt(params.sampleSize);
    firstFailureRate(idx) = revivalInfo.firstFailureRate;
    revivalCandidateRate(idx) = revivalInfo.candidateRate;
    revivalCandidateGivenFailure(idx) = revivalInfo.candidateGivenFailure;
    revivalObservedRate(idx) = revivalInfo.observedRate;
    revivalObservedGivenFailure(idx) = revivalInfo.observedGivenFailure;
    revivalObservedSE(idx) = revivalInfo.observedSE;
    meanRevivalDelay(idx) = revivalInfo.meanDelay;

    [availabilityTheory(idx), enumTotalProbability(idx)] = ...
        theoreticalAvailability(n, params.w, params);

    fprintf("%2d   %.6f   %.6f   %10.2f   %9.2f   %.6f   %.6f   %.3e\n", ...
        n, reliabilityMC(idx), reliabilitySE(idx), mttfMC(idx), ...
        mttfSE(idx), availabilityTheory(idx), revivalObservedGivenFailure(idx), ...
        abs(enumTotalProbability(idx) - 1));
end

if any(reliabilityMC < -eps | reliabilityMC > 1 + eps)
    error("Monte Carlo reliability is outside [0, 1].");
end
if any(availabilityTheory < -eps | availabilityTheory > 1 + eps)
    error("Theoretical availability is outside [0, 1].");
end
if max(abs(enumTotalProbability - 1)) > 1e-10
    error("The multinomial enumeration probability total is not close enough to 1.");
end

[bestReliability, bestReliabilityIndex] = max(reliabilityMC);
[bestMTTF, bestMTTFIndex] = max(mttfMC);
bestReliabilityN = params.nValues(bestReliabilityIndex);
bestMTTFN = params.nValues(bestMTTFIndex);

z95 = 1.96;
reliabilityCI95Low = max(0, reliabilityMC - z95 * reliabilitySE);
reliabilityCI95High = min(1, reliabilityMC + z95 * reliabilitySE);
mttfCI95Low = mttfMC - z95 * mttfSE;
mttfCI95High = mttfMC + z95 * mttfSE;
availabilityMinusReliability = availabilityTheory - reliabilityMC;
isBestReliability = params.nValues(:) == bestReliabilityN;
isBestMTTF = params.nValues(:) == bestMTTFN;

resultsTable = table(params.nValues(:), reliabilityMC, reliabilitySE, ...
    reliabilityCI95Low, reliabilityCI95High, mttfMC, mttfSE, ...
    mttfCI95Low, mttfCI95High, availabilityTheory, ...
    availabilityMinusReliability, firstFailureRate, revivalCandidateRate, ...
    revivalCandidateGivenFailure, revivalObservedRate, ...
    revivalObservedGivenFailure, revivalObservedSE, meanRevivalDelay, ...
    isBestReliability, isBestMTTF, abs(enumTotalProbability - 1), ...
    'VariableNames', {'n', 'R_MC_25000', 'R_MC_SE', ...
    'R_MC_CI95_low', 'R_MC_CI95_high', 'MTTF_MC', 'MTTF_MC_SE', ...
    'MTTF_CI95_low', 'MTTF_CI95_high', 'A_theory_25000', ...
    'A_minus_R', 'first_failure_rate', 'revival_candidate_rate', ...
    'revival_candidate_given_failure', 'revival_observed_rate', ...
    'revival_observed_given_failure', 'revival_observed_SE', ...
    'mean_revival_delay_hours', 'is_best_R', 'is_best_MTTF', ...
    'enum_total_abs_error'});

matPath = fullfile(resultDir, fileStem + "_results.mat");
csvPath = fullfile(resultDir, fileStem + "_results.csv");
summaryPath = fullfile(resultDir, fileStem + "_summary.txt");
reliabilityPlotPath = fullfile(figureDir, fileStem + "_reliability.png");
mttfPlotPath = fullfile(figureDir, fileStem + "_mttf.png");

save(matPath, "params", "resultsTable", "reliabilityMC", "reliabilitySE", ...
    "reliabilityCI95Low", "reliabilityCI95High", "mttfMC", "mttfSE", ...
    "mttfCI95Low", "mttfCI95High", "availabilityTheory", ...
    "availabilityMinusReliability", "firstFailureRate", ...
    "revivalCandidateRate", "revivalCandidateGivenFailure", ...
    "revivalObservedRate", "revivalObservedGivenFailure", ...
    "revivalObservedSE", "meanRevivalDelay", "enumTotalProbability", ...
    "bestReliabilityN", "bestReliability", "bestMTTFN", "bestMTTF", ...
    "isBestReliability", "isBestMTTF");
writetable(resultsTable, csvPath);

nValues = params.nValues(:);

fig = figure("Visible", "off", "Color", "w", "Position", [100, 100, 840, 520]);
ax = axes(fig);
hold(ax, "on");
reliabilityColor = [0.0000, 0.4470, 0.7410];
availabilityColor = [0.8500, 0.3250, 0.0980];
bestColor = [0.6350, 0.0780, 0.1840];
hReliability = errorbar(ax, nValues, reliabilityMC, z95 * reliabilitySE, ...
    "o-", "LineWidth", 1.9, "MarkerSize", 6, "CapSize", 7, ...
    "Color", reliabilityColor, "MarkerFaceColor", "w");
hAvailability = plot(ax, nValues, availabilityTheory, "s--", ...
    "LineWidth", 1.9, "MarkerSize", 6, "Color", availabilityColor, ...
    "MarkerFaceColor", "w");
hBest = plot(ax, bestReliabilityN, bestReliability, "p", ...
    "MarkerSize", 13, "MarkerFaceColor", bestColor, "MarkerEdgeColor", bestColor);
text(ax, bestReliabilityN + 0.25, bestReliability + 0.007, ...
    sprintf("best n = %d", bestReliabilityN), "Color", bestColor, ...
    "FontWeight", "bold", "FontSize", 10);
styleCase2Axes(ax, nValues);
ylim(ax, [0.32, 0.82]);
xlabel(ax, "number of nodes n", "FontSize", 11);
ylabel(ax, "probability at 25,000 hours", "FontSize", 11);
legend(ax, [hReliability, hAvailability, hBest], ...
    "Monte Carlo R(25,000) with 95% CI", ...
    "Theoretical availability A(25,000)", ...
    "Best reliability", "Location", "south", "Box", "on", "FontSize", 9);
exportgraphics(fig, reliabilityPlotPath, "Resolution", 300);
close(fig);

fig = figure("Visible", "off", "Color", "w", "Position", [100, 100, 840, 520]);
ax = axes(fig);
hold(ax, "on");
hMTTF = errorbar(ax, nValues, mttfMC / 1e4, z95 * mttfSE / 1e4, ...
    "o-", "LineWidth", 1.9, "MarkerSize", 6, "CapSize", 7, ...
    "Color", reliabilityColor, "MarkerFaceColor", "w");
hBestMTTF = plot(ax, bestMTTFN, bestMTTF / 1e4, "p", ...
    "MarkerSize", 13, "MarkerFaceColor", bestColor, "MarkerEdgeColor", bestColor);
text(ax, bestMTTFN + 0.25, bestMTTF / 1e4 + 0.07, ...
    sprintf("best n = %d", bestMTTFN), "Color", bestColor, ...
    "FontWeight", "bold", "FontSize", 10);
styleCase2Axes(ax, nValues);
ylim(ax, [2.35, 6.35]);
xlabel(ax, "number of nodes n", "FontSize", 11);
ylabel(ax, "MTTF (10^4 hours)", "FontSize", 11);
legend(ax, [hMTTF, hBestMTTF], "Monte Carlo MTTF with 95% CI", ...
    "Best MTTF", "Location", "northeast", "Box", "on", "FontSize", 9);
exportgraphics(fig, mttfPlotPath, "Resolution", 300);
close(fig);

runElapsed = toc(runTimer);
writeSummaryFile(summaryPath, params, bestReliabilityN, bestReliability, ...
    bestMTTFN, bestMTTF, resultsTable, runElapsed, matPath, csvPath, ...
    reliabilityPlotPath, mttfPlotPath);
save(matPath, "runElapsed", "summaryPath", "-append");

fprintf("\nBest R(25000): n = %d, R = %.6f\n", bestReliabilityN, bestReliability);
fprintf("Best MTTF:     n = %d, MTTF = %.2f hours\n", bestMTTFN, bestMTTF);
fprintf("Elapsed time: %.2f seconds\n", runElapsed);
fprintf("Saved:\n");
fprintf("  %s\n", matPath);
fprintf("  %s\n", csvPath);
fprintf("  %s\n", summaryPath);
fprintf("  %s\n", reliabilityPlotPath);
fprintf("  %s\n", mttfPlotPath);

function [lifetimes, revivalInfo] = simulateCase2Lifetimes(n, params)
    S = params.sampleSize;
    lifetimes = params.maxLife * ones(S, 1);
    firstRevivalCandidateTime = nan(S, 1);
    firstRevivalObservedTime = nan(S, 1);
    revivalCandidate = false(S, 1);
    revivalObserved = false(S, 1);
    revivalStream = RandStream("mt19937ar", "Seed", params.seed + 1009 + n);

    timeA = -params.meanA * log(rand(S, n));
    timeB = -params.meanB * log(rand(S, n));

    typeAUniform = rand(S, n);
    typeAFailure = ones(S, n, "uint8");
    typeAFailure(typeAUniform >= params.pAType(1) & ...
        typeAUniform < params.pAType(1) + params.pAType(2)) = 2;
    typeAFailure(typeAUniform >= params.pAType(1) + params.pAType(2)) = 3;

    typeBUniform = rand(S, n);
    typeBFailure = ones(S, n, "uint8");
    typeBFailure(typeBUniform >= params.pBType(1)) = 2;

    eventTimes = [timeA, timeB];
    [sortedTimes, eventOrder] = sort(eventTimes, 2);
    eventKind = [ones(1, n, "uint8"), 2 * ones(1, n, "uint8")];
    eventNode = [1:n, 1:n];

    stateA = zeros(S, n, "uint8");
    stateB = zeros(S, n, "uint8");
    nodeState = ones(S, n, "uint8");
    counts = zeros(S, 6);
    counts(:, 1) = n;
    alive = true(S, 1);

    for e = 1:(2 * n)
        rows = find(sortedTimes(:, e) <= params.maxLife);
        if isempty(rows)
            continue;
        end

        rows = rows(:);
        failedBeforeEvent = ~alive(rows);
        originalEvent = eventOrder(rows, e);
        originalEvent = originalEvent(:);
        nodes = eventNode(originalEvent);
        nodes = nodes(:);
        linearNodeIndex = sub2ind([S, n], rows, nodes);
        oldState = double(nodeState(linearNodeIndex));

        isAEvent = eventKind(originalEvent) == 1;
        isAEvent = isAEvent(:);
        if any(isAEvent)
            stateA(linearNodeIndex(isAEvent)) = typeAFailure(linearNodeIndex(isAEvent));
        end
        if any(~isAEvent)
            stateB(linearNodeIndex(~isAEvent)) = typeBFailure(linearNodeIndex(~isAEvent));
        end

        newState = mapNodeStateVector(stateA(linearNodeIndex), stateB(linearNodeIndex));
        changed = newState ~= oldState;
        if any(changed)
            changedRows = rows(changed);
            oldCountIndex = sub2ind([S, 6], changedRows, oldState(changed));
            newCountIndex = sub2ind([S, 6], changedRows, newState(changed));
            counts(oldCountIndex) = counts(oldCountIndex) - 1;
            counts(newCountIndex) = counts(newCountIndex) + 1;
            nodeState(linearNodeIndex(changed)) = uint8(newState(changed));
        end

        previouslyFailedRows = rows(failedBeforeEvent);
        if ~isempty(previouslyFailedRows)
            postFailureWorkProbability = systemWorkProbability( ...
                counts(previouslyFailedRows, :), params.k);

            newCandidate = postFailureWorkProbability > 0 & ...
                ~revivalCandidate(previouslyFailedRows);
            if any(newCandidate)
                candidateRows = previouslyFailedRows(newCandidate);
                revivalCandidate(candidateRows) = true;
                firstRevivalCandidateTime(candidateRows) = sortedTimes(candidateRows, e);
            end

            notYetObserved = ~revivalObserved(previouslyFailedRows);
            if any(notYetObserved)
                observationRows = previouslyFailedRows(notYetObserved);
                observationProbability = postFailureWorkProbability(notYetObserved);
                observedNow = observationProbability > 0 & ...
                    rand(revivalStream, numel(observationRows), 1) <= observationProbability;
                if any(observedNow)
                    revivedRows = observationRows(observedNow);
                    revivalObserved(revivedRows) = true;
                    firstRevivalObservedTime(revivedRows) = sortedTimes(revivedRows, e);
                end
            end
        end

        aliveRows = rows(~failedBeforeEvent);
        if isempty(aliveRows)
            continue;
        end

        workProbability = systemWorkProbability(counts(aliveRows, :), params.k);
        failedThisEvent = workProbability <= 0 | rand(numel(aliveRows), 1) > workProbability;
        if any(failedThisEvent)
            failedRows = aliveRows(failedThisEvent);
            lifetimes(failedRows) = sortedTimes(failedRows, e);
            alive(failedRows) = false;
        end
    end

    failedBeforeCap = lifetimes < params.maxLife;
    failedCount = sum(failedBeforeCap);
    if failedCount > 0
        observedGivenFailure = sum(revivalObserved) / failedCount;
        candidateGivenFailure = sum(revivalCandidate) / failedCount;
    else
        observedGivenFailure = 0;
        candidateGivenFailure = 0;
    end

    observedRate = mean(revivalObserved);
    delay = firstRevivalObservedTime - lifetimes;
    revivalInfo.firstFailureRate = mean(failedBeforeCap);
    revivalInfo.candidateRate = mean(revivalCandidate);
    revivalInfo.candidateGivenFailure = candidateGivenFailure;
    revivalInfo.observedRate = observedRate;
    revivalInfo.observedGivenFailure = observedGivenFailure;
    revivalInfo.observedSE = sqrt(observedRate * (1 - observedRate) / S);
    revivalInfo.meanDelay = mean(delay(revivalObserved), "omitnan");
end

function state = mapNodeStateVector(stateA, stateB)
    persistent stateMap
    if isempty(stateMap)
        % State codes: 1 PF, 2 MO, 3 SO, 4 FB, 5 DM, 6 DN.
        stateMap = uint8([ ...
            1, 2, 3; ...
            3, 4, 3; ...
            5, 2, 6; ...
            6, 6, 6]);
    end
    mapIndex = sub2ind([4, 3], double(stateA) + 1, double(stateB) + 1);
    state = double(stateMap(mapIndex));
end

function workProbability = systemWorkProbability(counts, k)
    qPF = counts(:, 1);
    qMO = counts(:, 2);
    qSO = counts(:, 3);
    qFB = counts(:, 4);
    qDM = counts(:, 5);

    workProbability = zeros(size(qPF));
    possible = qFB < 1 & qMO < 2 & qPF + qMO + qDM > 0 & ...
        qPF + qMO + qSO + qDM >= k;

    oneMO = possible & qMO == 1 & qPF + qSO >= k - 1;
    workProbability(oneMO) = 1;

    noMO = possible & qMO == 0;
    pfMasterWorks = noMO & qPF >= 1 & qPF + qSO >= k;
    dmOnlyMasterWorks = noMO & qPF == 0 & qDM >= 1 & qSO >= k - 1;
    criticalRandom = noMO & qPF >= 1 & qPF + qSO == k - 1 & qDM >= 1;

    workProbability(pfMasterWorks | dmOnlyMasterWorks) = 1;
    workProbability(criticalRandom) = qDM(criticalRandom) ./ ...
        (qDM(criticalRandom) + qPF(criticalRandom));
end

function [availability, totalProbability] = theoreticalAvailability(n, t, params)
    nodeProbability = nodeStateProbabilities(t, params);
    availability = 0;
    totalProbability = 0;

    for qPF = 0:n
        for qMO = 0:(n - qPF)
            for qSO = 0:(n - qPF - qMO)
                for qFB = 0:(n - qPF - qMO - qSO)
                    for qDM = 0:(n - qPF - qMO - qSO - qFB)
                        qDN = n - qPF - qMO - qSO - qFB - qDM;
                        counts = [qPF, qMO, qSO, qFB, qDM, qDN];

                        combinationProbability = multinomialProbability(counts, nodeProbability);
                        totalProbability = totalProbability + combinationProbability;
                        availability = availability + combinationProbability * ...
                            systemWorkProbability(counts, params.k);
                    end
                end
            end
        end
    end
end

function pNode = nodeStateProbabilities(t, params)
    pA0 = exp(-t / params.meanA);
    pA1 = params.pAType(1) * (1 - pA0);
    pA2 = params.pAType(2) * (1 - pA0);
    pA3 = params.pAType(3) * (1 - pA0);

    pB0 = exp(-t / params.meanB);
    pB1 = params.pBType(1) * (1 - pB0);
    pB2 = params.pBType(2) * (1 - pB0);

    pPF = pA0 * pB0;
    pMO = pA0 * pB1 + pA2 * pB1;
    pSO = pA0 * pB2 + pA1 * pB0 + pA1 * pB2;
    pFB = pA1 * pB1;
    pDM = pA2 * pB0;
    pDN = pA2 * pB2 + pA3 * (pB0 + pB1 + pB2);

    pNode = [pPF, pMO, pSO, pFB, pDM, pDN];
end

function probability = multinomialProbability(counts, categoryProbability)
    positive = counts > 0;
    if any(categoryProbability(positive) <= 0)
        probability = 0;
        return;
    end

    n = sum(counts);
    logProbability = gammaln(n + 1) - sum(gammaln(counts + 1)) + ...
        sum(counts(positive) .* log(categoryProbability(positive)));
    probability = exp(logProbability);
end

function styleCase2Axes(ax, nValues)
    grid(ax, "on");
    box(ax, "on");
    ax.FontName = "Arial";
    ax.FontSize = 10;
    ax.LineWidth = 0.9;
    ax.GridAlpha = 0.20;
    ax.TickDir = "out";
    ax.XTick = nValues(:)';
    ax.XLim = [min(nValues) - 0.25, max(nValues) + 0.25];
    ax.Units = "normalized";
    ax.Position = [0.10, 0.16, 0.86, 0.78];
end

function writeSummaryFile(summaryPath, params, bestReliabilityN, bestReliability, ...
    bestMTTFN, bestMTTF, resultsTable, runElapsed, matPath, csvPath, ...
    reliabilityPlotPath, mttfPlotPath)
    fid = fopen(summaryPath, "w", "n", "UTF-8");
    if fid < 0
        error("Failed to write summary file: %s", summaryPath);
    end

    cleanup = onCleanup(@() fclose(fid));
    bestRRow = resultsTable(resultsTable.is_best_R, :);
    bestMTTFRow = resultsTable(resultsTable.is_best_MTTF, :);
    [maxRevivalGivenFailure, maxRevivalIndex] = max(resultsTable.revival_observed_given_failure);
    maxRevivalRow = resultsTable(maxRevivalIndex, :);

    fprintf(fid, "Case 2 sonar reliability summary\n");
    fprintf(fid, "sample_size=%d\n", params.sampleSize);
    fprintf(fid, "seed=%d\n", params.seed);
    fprintf(fid, "k=%d\n", params.k);
    fprintf(fid, "w_hours=%.0f\n", params.w);
    fprintf(fid, "life_cap_hours=%.0f\n", params.maxLife);
    fprintf(fid, "mean_life_A_hours=%.0f\n", params.meanA);
    fprintf(fid, "mean_life_B_hours=%.0f\n", params.meanB);
    fprintf(fid, "best_R_n=%d\n", bestReliabilityN);
    fprintf(fid, "best_R_25000=%.6f\n", bestReliability);
    fprintf(fid, "best_R_CI95=[%.6f, %.6f]\n", ...
        bestRRow.R_MC_CI95_low, bestRRow.R_MC_CI95_high);
    fprintf(fid, "best_R_A_theory_25000=%.6f\n", bestRRow.A_theory_25000);
    fprintf(fid, "best_MTTF_n=%d\n", bestMTTFN);
    fprintf(fid, "best_MTTF_hours=%.2f\n", bestMTTF);
    fprintf(fid, "best_MTTF_CI95_hours=[%.2f, %.2f]\n", ...
        bestMTTFRow.MTTF_CI95_low, bestMTTFRow.MTTF_CI95_high);
    fprintf(fid, "R_at_best_MTTF=%.6f\n", bestMTTFRow.R_MC_25000);
    fprintf(fid, "revival_observed_given_failure_at_best_R=%.6f\n", ...
        bestRRow.revival_observed_given_failure);
    fprintf(fid, "revival_candidate_given_failure_at_best_R=%.6f\n", ...
        bestRRow.revival_candidate_given_failure);
    fprintf(fid, "mean_revival_delay_at_best_R_hours=%.2f\n", ...
        bestRRow.mean_revival_delay_hours);
    fprintf(fid, "max_revival_observed_given_failure_n=%d\n", maxRevivalRow.n);
    fprintf(fid, "max_revival_observed_given_failure=%.6f\n", maxRevivalGivenFailure);
    fprintf(fid, "elapsed_seconds=%.2f\n", runElapsed);
    fprintf(fid, "\nInterpretation:\n");
    fprintf(fid, "The Monte Carlo reliability R(25000) is maximized at n=%d.\n", bestReliabilityN);
    fprintf(fid, "The Monte Carlo MTTF point estimate is maximized at n=%d.\n", bestMTTFN);
    fprintf(fid, "The theoretical availability A(25000) is used only as an instantaneous-state approximation.\n");
    fprintf(fid, "The revival diagnostic continues each path after first failure and records whether a later state is workable again.\n");
    fprintf(fid, "Because revival after a first failure is possible, A(25000) can be higher than R(25000).\n");
    fprintf(fid, "\nOutput files:\n");
    fprintf(fid, "mat=%s\n", matPath);
    fprintf(fid, "csv=%s\n", csvPath);
    fprintf(fid, "reliability_plot=%s\n", reliabilityPlotPath);
    fprintf(fid, "mttf_plot=%s\n", mttfPlotPath);
end

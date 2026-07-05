function result = case1_evaluate_spline(dataset, scheme, options)
%case1_evaluate_spline Evaluate a scheme with MATLAB-style cubic spline.
%   This evaluator uses interp1(..., 'spline', 'extrap'), which corresponds
%   to MATLAB's not-a-knot cubic spline behavior. This is the convention
%   matching the baseline scheme P0=[4,19,38,54,79,88].

if nargin < 3
    options = struct();
end

if ~isfield(options, 'sampleIndices') || isempty(options.sampleIndices)
    sampleIndices = 1:dataset.sampleCount;
else
    sampleIndices = options.sampleIndices(:)';
end

if ~isfield(options, 'sampleOrder') || isempty(options.sampleOrder)
    orderedIndices = sampleIndices;
else
    orderedIndices = options.sampleOrder(:)';
    orderedIndices = orderedIndices(ismember(orderedIndices, sampleIndices));
    missingIndices = setdiff(sampleIndices, orderedIndices, 'stable');
    orderedIndices = [orderedIndices, missingIndices];
end

if ~isfield(options, 'costQ') || isempty(options.costQ)
    options.costQ = 70;
end

if ~isfield(options, 'needDetail') || isempty(options.needDetail)
    options.needDetail = false;
end

if ~isfield(options, 'enablePruning') || isempty(options.enablePruning)
    options.enablePruning = false;
end

if ~isfield(options, 'incumbentCost') || isempty(options.incumbentCost)
    options.incumbentCost = inf;
end

selectedIdx = normalize_scheme(dataset.pointCount, scheme);
selectedCount = numel(selectedIdx);
measurementCost = options.costQ * selectedCount;

temperatureGrid = dataset.temperatureGrid;
totalSampleCount = numel(orderedIndices);

sampleCost = zeros(totalSampleCount, 1);
processedCount = 0;
accumulatedCost = 0;
pruned = false;

if options.needDetail
    meanAbsError = zeros(totalSampleCount, 1);
    maxAbsError = zeros(totalSampleCount, 1);
    largeErrOver04 = zeros(totalSampleCount, 1);
    largeErrOver08 = zeros(totalSampleCount, 1);
    largeErrOver12 = zeros(totalSampleCount, 1);
    largeErrOver20 = zeros(totalSampleCount, 1);
end

for i = 1:totalSampleCount
    sampleId = orderedIndices(i);
    voltageRow = dataset.voltageMatrix(sampleId, :);
    voltageChosen = voltageRow(selectedIdx);
    temperatureChosen = temperatureGrid(selectedIdx);

    [voltageSorted, order] = sort(voltageChosen, 'ascend');
    temperatureSorted = temperatureChosen(order);
    [voltageUnique, uniqueIdx] = unique(voltageSorted, 'stable');
    temperatureUnique = temperatureSorted(uniqueIdx);

    if numel(voltageUnique) >= 2
        estimatedTemperature = interp1(voltageUnique, temperatureUnique, voltageRow, 'spline', 'extrap');
    else
        estimatedTemperature = zeros(size(voltageRow));
    end

    absError = abs(estimatedTemperature - temperatureGrid);
    errorCost = calc_error_cost(absError);
    oneSampleCost = sum(errorCost) + measurementCost;

    processedCount = processedCount + 1;
    sampleCost(processedCount) = oneSampleCost;
    accumulatedCost = accumulatedCost + oneSampleCost;

    if options.needDetail
        meanAbsError(processedCount) = mean(absError);
        maxAbsError(processedCount) = max(absError);
        largeErrOver04(processedCount) = sum(absError > 0.4);
        largeErrOver08(processedCount) = sum(absError > 0.8);
        largeErrOver12(processedCount) = sum(absError > 1.2);
        largeErrOver20(processedCount) = sum(absError > 2.0);
    end

    if options.enablePruning && isfinite(options.incumbentCost)
        remainingCount = totalSampleCount - processedCount;
        lowerBoundMeanCost = (accumulatedCost + remainingCount * measurementCost) / totalSampleCount;
        if lowerBoundMeanCost >= options.incumbentCost
            pruned = true;
            break;
        end
    end
end

sampleCost = sampleCost(1:processedCount);

if pruned
    remainingCount = totalSampleCount - processedCount;
    result.meanCost = (accumulatedCost + remainingCount * measurementCost) / totalSampleCount;
else
    result.meanCost = accumulatedCost / totalSampleCount;
end

result.selectedIdx = selectedIdx;
result.selectedTemperature = temperatureGrid(selectedIdx);
result.selectedCount = selectedCount;
result.sampleCost = sampleCost;
result.meanMeasurementCost = measurementCost;
result.meanErrorCost = result.meanCost - measurementCost;
result.processedSampleCount = processedCount;
result.totalSampleCount = totalSampleCount;
result.pruned = pruned;

if options.needDetail
    meanAbsError = meanAbsError(1:processedCount);
    maxAbsError = maxAbsError(1:processedCount);
    largeErrOver04 = largeErrOver04(1:processedCount);
    largeErrOver08 = largeErrOver08(1:processedCount);
    largeErrOver12 = largeErrOver12(1:processedCount);
    largeErrOver20 = largeErrOver20(1:processedCount);

    result.meanAbsError = mean(meanAbsError);
    result.maxAbsError = max(maxAbsError);
    result.avgCountOver04 = mean(largeErrOver04);
    result.avgCountOver08 = mean(largeErrOver08);
    result.avgCountOver12 = mean(largeErrOver12);
    result.avgCountOver20 = mean(largeErrOver20);
    result.totalCountOver20 = sum(largeErrOver20);
    result.sampleMeanAbsError = meanAbsError;
    result.sampleMaxAbsError = maxAbsError;
else
    result.meanAbsError = NaN;
    result.maxAbsError = NaN;
    result.avgCountOver04 = NaN;
    result.avgCountOver08 = NaN;
    result.avgCountOver12 = NaN;
    result.avgCountOver20 = NaN;
    result.totalCountOver20 = NaN;
    result.sampleMeanAbsError = [];
    result.sampleMaxAbsError = [];
end
end

function selectedIdx = normalize_scheme(pointCount, scheme)
scheme = scheme(:)';
if numel(scheme) == pointCount && all(ismember(unique(scheme), [0, 1]))
    selectedIdx = find(scheme > 0.5);
else
    selectedIdx = unique(round(scheme));
end
selectedIdx = selectedIdx(selectedIdx >= 1 & selectedIdx <= pointCount);
if isempty(selectedIdx)
    error('No calibration points are selected.');
end
end

function errorCost = calc_error_cost(absError)
errorCost = zeros(size(absError));
errorCost(absError > 0.4 & absError <= 0.8) = 2;
errorCost(absError > 0.8 & absError <= 1.2) = 15;
errorCost(absError > 1.2 & absError <= 1.6) = 30;
errorCost(absError > 1.6 & absError <= 2.0) = 45;
errorCost(absError > 2.0) = 80000;
end


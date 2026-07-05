function refineResult = case1_three_delete_three_add_refine(baseScheme, dataset, costQ, incumbentCost, method)
%case1_three_delete_three_add_refine Exact 3-delete-3-add search around one scheme.

if nargin < 4 || isempty(incumbentCost)
    incumbentCost = inf;
end
if nargin < 5 || isempty(method)
    method = 'spline';
end

baseScheme = sort(baseScheme(:)');
pointCount = numel(baseScheme);
removeTriplets = nchoosek(1:pointCount, 3);
allTriplets = nchoosek(1:dataset.pointCount, 3);
sampleOrder = 1:dataset.sampleCount;

bestScheme = baseScheme;
bestCost = incumbentCost;
bestSampleCost = [];
checkedCount = 0;

for rt = 1:size(removeTriplets, 1)
    reduced = baseScheme;
    reduced(removeTriplets(rt, :)) = [];

    validMask = ~ismember(allTriplets(:, 1), reduced) & ...
                ~ismember(allTriplets(:, 2), reduced) & ...
                ~ismember(allTriplets(:, 3), reduced);
    addTriplets = allTriplets(validMask, :);

    for at = 1:size(addTriplets, 1)
        checkedCount = checkedCount + 1;
        if mod(checkedCount, 50000) == 0
            fprintf('[3del3add] checked %d, current best = %.3f\n', checkedCount, bestCost);
        end

        candidate = sort([reduced, addTriplets(at, :)]);

        evalOptions.sampleIndices = 1:dataset.sampleCount;
        evalOptions.sampleOrder = sampleOrder;
        evalOptions.costQ = costQ;
        evalOptions.needDetail = false;
        evalOptions.enablePruning = true;
        evalOptions.incumbentCost = bestCost;
        evalOptions.method = method;

        evalResult = case1_evaluate_scheme(dataset, candidate, evalOptions);

        if evalResult.meanCost < bestCost
            bestCost = evalResult.meanCost;
            bestScheme = candidate;
            bestSampleCost = evalResult.sampleCost;
            fprintf('[3del3add] checked %d, improved to %.3f, idx = %s\n', ...
                checkedCount, bestCost, mat2str(bestScheme));
        end
    end
end

refineResult.bestScheme = bestScheme;
refineResult.bestCost = bestCost;
refineResult.sampleCost = bestSampleCost;
refineResult.checkedCount = checkedCount;
end


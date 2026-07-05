function refineResult = case1_two_delete_two_add_refine(baseScheme, dataset, costQ, incumbentCost, method)
%case1_two_delete_two_add_refine Exact 2-delete-2-add search around one scheme.

if nargin < 4 || isempty(incumbentCost)
    incumbentCost = inf;
end
if nargin < 5 || isempty(method)
    method = 'spline';
end

baseScheme = sort(baseScheme(:)');
pointCount = numel(baseScheme);
removePairs = nchoosek(1:pointCount, 2);
sampleOrder = 1:dataset.sampleCount;

bestScheme = baseScheme;
bestCost = incumbentCost;
bestSampleCost = [];
checkedCount = 0;

for rp = 1:size(removePairs, 1)
    reduced = baseScheme;
    reduced(removePairs(rp, :)) = [];

    addPool = setdiff(1:dataset.pointCount, reduced);
    addPairs = nchoosek(addPool, 2);

    for ap = 1:size(addPairs, 1)
        checkedCount = checkedCount + 1;
        candidate = sort([reduced, addPairs(ap, :)]);

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
            fprintf('[2del2add] checked %d, improved to %.3f, idx = %s\n', ...
                checkedCount, bestCost, mat2str(bestScheme));
        end
    end
end

refineResult.bestScheme = bestScheme;
refineResult.bestCost = bestCost;
refineResult.sampleCost = bestSampleCost;
refineResult.checkedCount = checkedCount;
end


function refineResult = case1_refine_scheme_local(dataset, initialScheme, sampleIndices, sampleOrder, options)
%case1_refine_scheme_local Lightweight local search using v7 objective.

if nargin < 5
    options = struct();
end

if ~isfield(options, 'window') || isempty(options.window)
    options.window = 2;
end
if ~isfield(options, 'maxPasses') || isempty(options.maxPasses)
    options.maxPasses = 1;
end
if ~isfield(options, 'costQ') || isempty(options.costQ)
    options.costQ = 70;
end
if ~isfield(options, 'enablePruning') || isempty(options.enablePruning)
    options.enablePruning = true;
end
if ~isfield(options, 'method') || isempty(options.method)
    options.method = 'spline';
end

geneCount = dataset.pointCount;
pointCount = numel(initialScheme);
cache = containers.Map('KeyType', 'char', 'ValueType', 'double');

currentScheme = repair_scheme(initialScheme, geneCount, pointCount);
currentCost = evaluate_scheme_cost(currentScheme, inf);

for passId = 1:options.maxPasses
    improved = false;
    bestScheme = currentScheme;
    bestCost = currentCost;

    for geneId = 1:pointCount
        lowerBound = 1;
        upperBound = geneCount;
        if geneId > 1
            lowerBound = currentScheme(geneId - 1) + 1;
        end
        if geneId < pointCount
            upperBound = currentScheme(geneId + 1) - 1;
        end

        candidateRange = max(lowerBound, currentScheme(geneId) - options.window): ...
            min(upperBound, currentScheme(geneId) + options.window);
        candidateRange(candidateRange == currentScheme(geneId)) = [];

        for newIndex = candidateRange
            candidate = currentScheme;
            candidate(geneId) = newIndex;
            candidate = repair_scheme(candidate, geneCount, pointCount);
            candidateCost = evaluate_scheme_cost(candidate, bestCost);
            if candidateCost + 1e-9 < bestCost
                bestScheme = candidate;
                bestCost = candidateCost;
                improved = true;
            end
        end
    end

    if improved
        currentScheme = bestScheme;
        currentCost = bestCost;
    else
        break;
    end
end

refineResult.bestScheme = currentScheme;
refineResult.bestCost = currentCost;

    function costValue = evaluate_scheme_cost(scheme, incumbentCost)
        key = sprintf('%03d_', scheme);
        if isKey(cache, key)
            costValue = cache(key);
            return;
        end
        evalOptions.sampleIndices = sampleIndices;
        evalOptions.sampleOrder = sampleOrder;
        evalOptions.costQ = options.costQ;
        evalOptions.needDetail = false;
        evalOptions.enablePruning = options.enablePruning;
        evalOptions.incumbentCost = incumbentCost;
        evalOptions.method = options.method;
        evalResult = case1_evaluate_scheme(dataset, scheme, evalOptions);
        costValue = evalResult.meanCost;
        cache(key) = costValue;
    end
end

function scheme = repair_scheme(rawScheme, geneCount, pointCount)
rawScheme = unique(round(rawScheme(:)'));
rawScheme(rawScheme < 1 | rawScheme > geneCount) = [];
scheme = sort(rawScheme);
while numel(scheme) < pointCount
    pool = setdiff(1:geneCount, scheme);
    if isempty(pool)
        break;
    end
    scheme(end + 1) = pool(randi(numel(pool))); %#ok<AGROW>
    scheme = sort(unique(scheme));
end
if numel(scheme) > pointCount
    scheme = scheme(1:pointCount);
end
end


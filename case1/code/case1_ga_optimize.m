function optimizerResult = case1_ga_optimize(dataset, cfg, resultDir)
%case1_ga_optimize Variable-point-count memetic GA for Case1.

search = cfg.search;
evalCfg = cfg.eval;
localCfg = cfg.localSearch;
randomCfg = cfg.random;
outputCfg = cfg.output;

geneCount = dataset.pointCount;
countValues = search.minPointCount:search.maxPointCount;
priority = build_priority_profile_SUBMIT(dataset);

restartSeeds = randomCfg.baseSeed + (0:randomCfg.restartCount - 1);
totalIterations = randomCfg.restartCount * search.maxGenerations;
fullSampleIndices = 1:dataset.sampleCount;
baseSearchSampleIndices = build_search_sample_indices_SUBMIT(dataset.sampleCount, evalCfg.searchSampleCount, ...
    evalCfg.searchRandomSampleCount, randomCfg.baseSeed);

exactCache = containers.Map('KeyType', 'char', 'ValueType', 'any');
archive = struct('scheme', {}, 'cost', {}, 'count', {}, 'selectedIdx', {}, 'temperature', {});
bestByCount = repmat(empty_record_SUBMIT(), numel(countValues), 1);
for i = 1:numel(countValues)
    bestByCount(i).count = countValues(i);
end
globalBest = empty_record_SUBMIT();
hardSamplePool = [];

historyRows = zeros(totalIterations, 12);
historyCursor = 0;
startTimer = tic;

for restartId = 1:numel(restartSeeds)
    rng(restartSeeds(restartId), 'twister');
    islands = initialize_islands_SUBMIT(dataset, search, priority);
    stagnationCount = 0;
    lastBestCost = globalBest.cost;

    for gen = 1:search.maxGenerations
        iteration = (restartId - 1) * search.maxGenerations + gen;
        sampleOrder = build_priority_order_SUBMIT(fullSampleIndices, hardSamplePool);
        searchSampleIndices = build_iteration_sample_indices_SUBMIT( ...
            baseSearchSampleIndices, hardSamplePool, fullSampleIndices, evalCfg.searchSampleCount);
        searchSampleOrder = build_priority_order_SUBMIT(searchSampleIndices, hardSamplePool);
        adaptive = adaptive_params_SUBMIT(search, gen, search.maxGenerations, stagnationCount, islands);

        islandSummaries = repmat(struct( ...
            'population', [], ...
            'searchCost', [], ...
            'exactCandidates', [], ...
            'bestSearchCost', inf, ...
            'bestSearchCount', NaN, ...
            'meanCount', NaN, ...
            'diversity', NaN, ...
            'prunedRatio', NaN), search.islandCount, 1);

        exactCandidatePool = zeros(search.islandCount * max(search.exactEliteCount, 1) + ...
            max(1, search.archiveInjectCount) + 2, geneCount);
        exactCursor = 0;

        generationBestSearchCost = inf;
        generationBestSearchCount = NaN;
        generationSearchCostSum = 0;
        generationSearchCostCount = 0;
        meanIslandDiversity = 0;
        meanCountValue = 0;
        meanPrunedRatio = 0;

        for islandId = 1:search.islandCount
            population = islands{islandId};
            [population, searchCost, exactCandidates, meta, exactCache] = evaluate_island_population_SUBMIT( ...
                population, dataset, searchSampleIndices, searchSampleOrder, globalBest.cost, cfg, priority, exactCache);

            islandSummaries(islandId).population = population;
            islandSummaries(islandId).searchCost = searchCost;
            islandSummaries(islandId).exactCandidates = exactCandidates;
            islandSummaries(islandId).bestSearchCost = searchCost(1);
            islandSummaries(islandId).bestSearchCount = nnz(population(1, :));
            islandSummaries(islandId).meanCount = mean(sum(population, 2));
            islandSummaries(islandId).diversity = meta.diversity;
            islandSummaries(islandId).prunedRatio = meta.prunedRatio;

            exactCount = min(size(exactCandidates, 1), search.exactEliteCount);
            exactCandidatePool(exactCursor + 1:exactCursor + exactCount, :) = exactCandidates(1:exactCount, :);
            exactCursor = exactCursor + exactCount;

            generationBestSearchCost = min(generationBestSearchCost, searchCost(1));
            generationSearchCostSum = generationSearchCostSum + sum(searchCost);
            generationSearchCostCount = generationSearchCostCount + numel(searchCost);
            if islandSummaries(islandId).bestSearchCost <= generationBestSearchCost + 1e-9
                generationBestSearchCount = islandSummaries(islandId).bestSearchCount;
            end
            meanIslandDiversity = meanIslandDiversity + islandSummaries(islandId).diversity;
            meanCountValue = meanCountValue + islandSummaries(islandId).meanCount;
            meanPrunedRatio = meanPrunedRatio + islandSummaries(islandId).prunedRatio;
        end

        if ~isempty(globalBest.scheme)
            exactCursor = exactCursor + 1;
            exactCandidatePool(exactCursor, :) = globalBest.scheme;
        end
        for archiveId = 1:min(numel(archive), search.archiveInjectCount)
            exactCursor = exactCursor + 1;
            exactCandidatePool(exactCursor, :) = archive(archiveId).scheme;
        end
        exactCandidatePool = exactCandidatePool(1:exactCursor, :);
        exactCandidatePool = unique(exactCandidatePool, 'rows', 'stable');

        [generationBestExact, globalBest, bestByCount, archive, exactCache] = exact_update_candidates_SUBMIT( ...
            exactCandidatePool, dataset, sampleOrder, cfg, priority, globalBest, bestByCount, archive, exactCache);

        if localCfg.enable && ~isempty(globalBest.scheme)
            shouldRefine = gen == 1 || mod(gen, localCfg.every) == 0 || stagnationCount >= search.stagnationWindow;
            if shouldRefine
                localOptions.maxPasses = localCfg.maxPasses;
                localOptions.moveWindow = localCfg.moveWindow;
                localOptions.addCandidateCount = localCfg.addCandidateCount;
                localOptions.removeCandidateCount = localCfg.removeCandidateCount;
                localOptions.swapAddCandidateCount = localCfg.swapAddCandidateCount;
                localOptions.swapRemoveCandidateCount = localCfg.swapRemoveCandidateCount;
                [localRecord, exactCache] = refine_scheme_variable_SUBMIT( ...
                    dataset, globalBest, sampleOrder, cfg, priority, exactCache, localOptions);
                if localRecord.cost < globalBest.cost
                    [globalBest, bestByCount, archive] = assimilate_exact_records_SUBMIT( ...
                        localRecord, globalBest, bestByCount, archive, cfg);
                end
            end
        end

        if ~isempty(globalBest.sampleCost)
            [~, hardOrder] = sort(globalBest.sampleCost, 'descend');
            hardSamplePool = hardOrder(1:min(evalCfg.hardSampleCount, numel(hardOrder)));
        end

        if globalBest.cost + 1e-9 < lastBestCost
            stagnationCount = 0;
            lastBestCost = globalBest.cost;
        else
            stagnationCount = stagnationCount + 1;
        end

        islands = reproduce_islands_SUBMIT(islandSummaries, dataset, cfg, priority, archive, adaptive);
        if mod(gen, search.migrationEvery) == 0
            islands = migrate_islands_SUBMIT(islands, islandSummaries, search.migrationCount);
        end

        historyCursor = historyCursor + 1;
        elapsedSeconds = toc(startTimer);
        generationMeanSearchCost = generationSearchCostSum / max(generationSearchCostCount, 1);
        historyRows(historyCursor, :) = [ ...
            restartId, gen, iteration, generationBestSearchCost, generationMeanSearchCost, generationBestExact, globalBest.cost, ...
            globalBest.count, meanIslandDiversity / search.islandCount, meanCountValue / search.islandCount, ...
            meanPrunedRatio / search.islandCount, elapsedSeconds];

        if outputCfg.consoleEvery > 0 && ...
                (gen == 1 || mod(gen, outputCfg.consoleEvery) == 0 || gen == search.maxGenerations)
            fprintf(['[SUBMIT][R %d/%d][G %d/%d] search best = %.3f (k=%d), search mean = %.3f, ', ...
                'exact best = %.3f (k=%d), islandDiv = %.3f, meanCount = %.2f, pruned = %.2f, elapsed = %.1fs\n'], ...
                restartId, numel(restartSeeds), gen, search.maxGenerations, ...
                generationBestSearchCost, generationBestSearchCount, generationMeanSearchCost, ...
                globalBest.cost, globalBest.count, ...
                meanIslandDiversity / search.islandCount, meanCountValue / search.islandCount, ...
                meanPrunedRatio / search.islandCount, elapsedSeconds);
        end

        if outputCfg.checkpointEvery > 0 && mod(iteration, outputCfg.checkpointEvery) == 0
            checkpoint.globalBest = globalBest;
            checkpoint.bestByCount = bestByCount;
            checkpoint.archive = archive;
            checkpoint.historyRows = historyRows(1:historyCursor, :);
            save(fullfile(resultDir, 'latest_checkpoint.mat'), 'checkpoint');
        end
    end
end

if localCfg.enable && ~isempty(globalBest.scheme)
    sampleOrder = build_priority_order_SUBMIT(fullSampleIndices, hardSamplePool);
    localOptions.maxPasses = localCfg.finalPasses;
    localOptions.moveWindow = localCfg.finalMoveWindow;
    localOptions.addCandidateCount = localCfg.finalAddCandidateCount;
    localOptions.removeCandidateCount = localCfg.finalRemoveCandidateCount;
    localOptions.swapAddCandidateCount = localCfg.finalAddCandidateCount;
    localOptions.swapRemoveCandidateCount = localCfg.finalRemoveCandidateCount;
    [localRecord, exactCache] = refine_scheme_variable_SUBMIT( ...
        dataset, globalBest, sampleOrder, cfg, priority, exactCache, localOptions);
    if localRecord.cost < globalBest.cost
        [globalBest, bestByCount, archive] = assimilate_exact_records_SUBMIT( ...
            localRecord, globalBest, bestByCount, archive, cfg);
    end
end

transitionRecords = repmat(empty_record_SUBMIT(), 0, 1);
if cfg.finalRefine.enableCountTransition && ~isempty(globalBest.scheme)
    sampleOrder = build_priority_order_SUBMIT(fullSampleIndices, hardSamplePool);
    [transitionRecords, exactCache] = count_transition_refine_SUBMIT( ...
        dataset, globalBest, bestByCount, archive, sampleOrder, cfg, priority, exactCache);
    if ~isempty(transitionRecords)
        [globalBest, bestByCount, archive] = assimilate_exact_records_SUBMIT( ...
            transitionRecords, globalBest, bestByCount, archive, cfg);
    end
end

if cfg.finalRefine.enableTwoDeleteTwoAdd
    sampleOrder = 1:dataset.sampleCount;
    [twoDeleteRecords, exactCache] = iterative_two_delete_refine_SUBMIT( ...
        dataset, globalBest, bestByCount, archive, transitionRecords, sampleOrder, cfg, priority, exactCache);
    if ~isempty(twoDeleteRecords)
        [globalBest, bestByCount, archive] = assimilate_exact_records_SUBMIT( ...
            twoDeleteRecords, globalBest, bestByCount, archive, cfg);
    end
end

if cfg.finalRefine.enableThreeDeleteThreeAdd && globalBest.count <= 6
    exactRefine3 = case1_three_delete_three_add_refine(globalBest.selectedIdx, dataset, cfg.eval.costQ, globalBest.cost, cfg.eval.method);
    if exactRefine3.bestCost < globalBest.cost
        sampleOrder = 1:dataset.sampleCount;
        [exactRecord, exactCache] = exact_evaluate_scheme_SUBMIT( ...
            scheme_from_indices_SUBMIT(exactRefine3.bestScheme, geneCount), dataset, sampleOrder, cfg, priority, exactCache);
        [globalBest, bestByCount, archive] = assimilate_exact_records_SUBMIT( ...
            exactRecord, globalBest, bestByCount, archive, cfg);
    end
end

optimizerResult.bestScheme = globalBest.scheme;
optimizerResult.bestCost = globalBest.cost;
optimizerResult.bestCount = globalBest.count;
optimizerResult.bestSelectedIdx = globalBest.selectedIdx;
optimizerResult.bestTemperature = globalBest.temperature;
optimizerResult.bestByCount = bestByCount;
optimizerResult.archive = archive;
optimizerResult.historyTable = array2table(historyRows(1:historyCursor, :), 'VariableNames', { ...
    'Restart', 'Generation', 'Iteration', 'GenerationBestSearchCost', 'GenerationMeanSearchCost', 'GenerationBestExactCost', ...
    'GlobalBestExactCost', 'GlobalBestCount', 'IslandDiversity', 'MeanPopulationCount', ...
    'PrunedRatio', 'ElapsedSeconds'});
end

function [records, exactCache] = iterative_two_delete_refine_SUBMIT( ...
        dataset, globalBest, bestByCount, archive, transitionRecords, sampleOrder, cfg, priority, exactCache)
records = repmat(empty_record_SUBMIT(), 0, 1);
sources = collect_two_delete_sources_SUBMIT(globalBest, bestByCount, archive, transitionRecords, cfg);

for sourceId = 1:numel(sources)
    source = sources(sourceId);
    if isempty(source.scheme) || source.count ~= 6 || source.cost > cfg.finalRefine.twoDeleteMaxStartCost
        continue;
    end

    currentScheme = source.selectedIdx;
    currentCost = source.cost;
    fprintf('[SUBMIT 2del2add] source %d: cost=%.3f, T=%s\n', ...
        sourceId, currentCost, vector_to_string_SUBMIT(source.temperature));

    for passId = 1:cfg.finalRefine.twoDeletePasses
        refineResult = case1_two_delete_two_add_refine( ...
            currentScheme, dataset, cfg.eval.costQ, currentCost, cfg.eval.method);
        if refineResult.bestCost >= currentCost - 1e-9
            break;
        end

        currentScheme = refineResult.bestScheme;
        currentCost = refineResult.bestCost;
        [record, exactCache] = exact_evaluate_scheme_SUBMIT( ...
            scheme_from_indices_SUBMIT(currentScheme, dataset.pointCount), dataset, sampleOrder, cfg, priority, exactCache);
        records(end + 1, 1) = record; %#ok<AGROW>

        fprintf('[SUBMIT 2del2add] pass %d -> cost=%.3f, T=%s\n', ...
            passId, record.cost, vector_to_string_SUBMIT(record.temperature));
    end
end
end

function sources = collect_two_delete_sources_SUBMIT(globalBest, bestByCount, archive, transitionRecords, cfg)
rawSources = repmat(empty_record_SUBMIT(), 0, 1);
if ~isempty(globalBest.scheme)
    rawSources(end + 1, 1) = globalBest; %#ok<AGROW>
end

for i = 1:numel(bestByCount)
    if ~isempty(bestByCount(i).scheme)
        rawSources(end + 1, 1) = bestByCount(i); %#ok<AGROW>
    end
end

for i = 1:numel(archive)
    rawSources(end + 1, 1) = archive(i); %#ok<AGROW>
end

for i = 1:numel(transitionRecords)
    if ~isempty(transitionRecords(i).scheme)
        rawSources(end + 1, 1) = transitionRecords(i); %#ok<AGROW>
    end
end

rawSources = rawSources(arrayfun(@(x) ~isempty(x.scheme) && x.count == 6 && ...
    x.cost <= cfg.finalRefine.twoDeleteMaxStartCost, rawSources));
if isempty(rawSources)
    sources = rawSources;
    return;
end

[~, order] = sort([rawSources.cost], 'ascend');
rawSources = rawSources(order);

sources = repmat(empty_record_SUBMIT(), 0, 1);
for i = 1:numel(rawSources)
    candidate = rawSources(i);
    duplicate = false;
    for j = 1:numel(sources)
        if archive_distance_SUBMIT(candidate.scheme, sources(j).scheme) < 4
            duplicate = true;
            break;
        end
    end
    if duplicate
        continue;
    end

    sources(end + 1, 1) = candidate; %#ok<AGROW>
    if numel(sources) >= cfg.finalRefine.twoDeleteSourceCount
        break;
    end
end
end

function [records, exactCache] = count_transition_refine_SUBMIT( ...
        dataset, globalBest, bestByCount, archive, sampleOrder, cfg, priority, exactCache)
records = repmat(empty_record_SUBMIT(), 0, 1);
sourceRecords = collect_transition_sources_SUBMIT(globalBest, bestByCount, archive, cfg.finalRefine.countTransitionTopCount);

if isempty(sourceRecords)
    return;
end

seenKeys = containers.Map('KeyType', 'char', 'ValueType', 'logical');
usedSourceCount = 0;

for sourceId = 1:numel(sourceRecords)
    source = sourceRecords(sourceId);
    if isempty(source.scheme) || source.count <= 6 || ...
            source.count > cfg.finalRefine.countTransitionMaxSourceCount
        continue;
    end

    sourceKey = scheme_key_SUBMIT(source.scheme);
    if isKey(seenKeys, sourceKey)
        continue;
    end
    seenKeys(sourceKey) = true;
    usedSourceCount = usedSourceCount + 1;

    selectedIdx = source.selectedIdx;
    if isempty(selectedIdx)
        selectedIdx = indices_from_scheme_SUBMIT(source.scheme);
    end

    fprintf('[SUBMIT transition] source %d: k=%d cost=%.3f, T=%s\n', ...
        usedSourceCount, source.count, source.cost, vector_to_string_SUBMIT(source.temperature));

    for dropId = 1:numel(selectedIdx)
        seed = selectedIdx;
        seed(dropId) = [];
        if numel(seed) < cfg.search.minPointCount
            continue;
        end

        localOptions.window = cfg.finalRefine.countTransitionWindow;
        localOptions.maxPasses = cfg.finalRefine.countTransitionPasses;
        localOptions.costQ = cfg.eval.costQ;
        localOptions.enablePruning = true;
        localOptions.method = cfg.eval.method;

        refineResult = case1_refine_scheme_local( ...
            dataset, seed, 1:dataset.sampleCount, sampleOrder, localOptions);
        candidateScheme = scheme_from_indices_SUBMIT(refineResult.bestScheme, dataset.pointCount);
        [record, exactCache] = exact_evaluate_scheme_SUBMIT( ...
            candidateScheme, dataset, sampleOrder, cfg, priority, exactCache);
        records(end + 1, 1) = record; %#ok<AGROW>

        fprintf('[SUBMIT transition] drop %d -> k=%d cost=%.3f, T=%s\n', ...
            dropId, record.count, record.cost, vector_to_string_SUBMIT(record.temperature));
    end

    if usedSourceCount >= cfg.finalRefine.countTransitionTopCount
        break;
    end
end
end

function sourceRecords = collect_transition_sources_SUBMIT(globalBest, bestByCount, archive, topCount)
sourceRecords = repmat(empty_record_SUBMIT(), 0, 1);
if ~isempty(globalBest.scheme)
    sourceRecords(end + 1, 1) = globalBest; %#ok<AGROW>
end

for i = 1:numel(bestByCount)
    if ~isempty(bestByCount(i).scheme)
        sourceRecords(end + 1, 1) = bestByCount(i); %#ok<AGROW>
    end
end

archiveCount = min(numel(archive), max(1, topCount));
for i = 1:archiveCount
    sourceRecords(end + 1, 1) = archive(i); %#ok<AGROW>
end

if isempty(sourceRecords)
    return;
end

[~, order] = sort([sourceRecords.cost], 'ascend');
sourceRecords = sourceRecords(order);
end

function priority = build_priority_profile_SUBMIT(dataset)
meanVoltage = mean(dataset.voltageMatrix, 1);
slope = [abs(diff(meanVoltage)), 0];
curvature = [0, abs(diff(meanVoltage, 2)), 0];
sampleStd = std(dataset.voltageMatrix, 0, 1);
score = normalize01_SUBMIT(curvature) + 0.70 * normalize01_SUBMIT(sampleStd) + ...
    0.35 * normalize01_SUBMIT(slope) + 0.05;
score = score(:)';

[~, order] = sort(score, 'descend');
priority.score = score;
priority.order = order;
priority.meanVoltage = meanVoltage;
priority.curvature = curvature;
priority.sampleStd = sampleStd;
end

function islands = initialize_islands_SUBMIT(dataset, search, priority)
basePopulation = initialize_population_SUBMIT(dataset, search, priority);
islands = cell(search.islandCount, 1);
for islandId = 1:search.islandCount
    shifted = basePopulation;
    for rowId = 1:size(shifted, 1)
        shifted(rowId, :) = mutate_binary_scheme_SUBMIT( ...
            shifted(rowId, :), dataset.pointCount, search, priority, 0.12, 2);
    end
    islands{islandId} = shifted;
end
end

function population = initialize_population_SUBMIT(dataset, search, priority)
populationSize = search.populationSizePerIsland;
seedSet = build_neutral_seed_schemes_SUBMIT(dataset, search, priority);
seedSet = unique(seedSet, 'rows', 'stable');

seedKeepCount = min(size(seedSet, 1), max(1, ceil(0.55 * populationSize)));
if size(seedSet, 1) > seedKeepCount
    seedRows = unique(round(linspace(1, size(seedSet, 1), seedKeepCount)), 'stable');
    seedSet = seedSet(seedRows, :);
end

population = seedSet;
if size(population, 1) < populationSize
    extra = random_binary_population_SUBMIT(populationSize - size(population, 1), dataset.pointCount, search, priority);
    population = [population; extra]; %#ok<AGROW>
end
population = unique(population, 'rows', 'stable');

while size(population, 1) < populationSize
    newcomer = random_binary_population_SUBMIT(1, dataset.pointCount, search, priority);
    population = unique([population; newcomer], 'rows', 'stable'); %#ok<AGROW>
end

population = population(1:populationSize, :);
end

function seedSet = build_neutral_seed_schemes_SUBMIT(dataset, search, priority)
geneCount = dataset.pointCount;
countValues = search.minPointCount:search.maxPointCount;
seedSet = zeros(0, geneCount);

for pointCount = countValues
    idx1 = round(linspace(1, geneCount, pointCount));
    idx2 = round(((1:pointCount) - 0.5) * (geneCount / pointCount));
    idx2 = min(geneCount, max(1, idx2));
    idx3 = sort(priority.order(1:pointCount));
    idx4 = build_segment_seed_SUBMIT(priority.score, pointCount);
    idx5 = sort(unique(round(linspace(2, geneCount - 1, pointCount))));

    rawSeeds = {
        idx1;
        idx2;
        idx3;
        idx4;
        idx5};

    for seedId = 1:numel(rawSeeds)
        repaired = repair_binary_scheme_SUBMIT( ...
            scheme_from_indices_SUBMIT(rawSeeds{seedId}, geneCount), search, priority);
        seedSet(end + 1, :) = repaired; %#ok<AGROW>
    end
end

for jitterId = 1:search.seedJitterCount
    sourceId = randi(size(seedSet, 1));
    scheme = seedSet(sourceId, :);
    selectedIdx = find(scheme);
    if isempty(selectedIdx)
        continue;
    end
    replaceCount = min(numel(selectedIdx), randi([1, 2]));
    targetIdx = selectedIdx(randperm(numel(selectedIdx), replaceCount));
    newScheme = scheme;
    for pickId = 1:numel(targetIdx)
        oldPos = targetIdx(pickId);
        candidate = oldPos + randi([-search.seedJitterRadius, search.seedJitterRadius]);
        candidate = min(geneCount, max(1, candidate));
        newScheme(oldPos) = 0;
        newScheme(candidate) = 1;
    end
    newScheme = repair_binary_scheme_SUBMIT(newScheme, search, priority);
    seedSet(end + 1, :) = newScheme; %#ok<AGROW>
end
end

function idx = build_segment_seed_SUBMIT(score, pointCount)
geneCount = numel(score);
edges = round(linspace(1, geneCount + 1, pointCount + 1));
idx = zeros(1, pointCount);

for segId = 1:pointCount
    left = edges(segId);
    right = edges(segId + 1) - 1;
    right = max(left, min(geneCount, right));
    segment = left:right;
    [~, localId] = max(score(segment));
    idx(segId) = segment(localId);
end

idx = sort(unique(idx));
end

function [population, searchCost, exactCandidates, meta, exactCache] = evaluate_island_population_SUBMIT( ...
        population, dataset, sampleIndices, sampleOrder, incumbentCost, cfg, priority, exactCache)
population = unique(population, 'rows', 'stable');
memberCount = size(population, 1);

searchCost = inf(memberCount, 1);
processedCount = zeros(memberCount, 1);
prunedFlag = false(memberCount, 1);
exactMask = false(memberCount, 1);

for rowId = 1:memberCount
    [costValue, exactFlag, processedSamples, wasPruned, exactCache] = search_evaluate_scheme_SUBMIT( ...
        population(rowId, :), dataset, sampleIndices, sampleOrder, incumbentCost, cfg, priority, exactCache);
    searchCost(rowId) = costValue;
    processedCount(rowId) = processedSamples;
    prunedFlag(rowId) = wasPruned;
    exactMask(rowId) = exactFlag;

    if exactFlag && costValue < incumbentCost
        incumbentCost = costValue;
    end
end

[searchCost, order] = sort(searchCost, 'ascend');
population = population(order, :);
processedCount = processedCount(order);
prunedFlag = prunedFlag(order);
exactMask = exactMask(order);

exactCandidates = population;
meta.prunedRatio = mean(prunedFlag);
meta.diversity = compute_population_diversity_SUBMIT(population);
meta.meanProcessedSamples = mean(processedCount);
meta.exactRatio = mean(exactMask);
end

function [costValue, exactFlag, processedSamples, wasPruned, exactCache] = search_evaluate_scheme_SUBMIT( ...
        scheme, dataset, sampleIndices, sampleOrder, incumbentCost, cfg, priority, exactCache)
scheme = repair_binary_scheme_SUBMIT(scheme, cfg.search, priority);
key = scheme_key_SUBMIT(scheme);

if isKey(exactCache, key)
    record = exactCache(key);
    costValue = record.cost;
    exactFlag = true;
    processedSamples = dataset.sampleCount;
    wasPruned = false;
    return;
end

evalOptions.sampleIndices = sampleIndices;
evalOptions.sampleOrder = sampleOrder;
evalOptions.costQ = cfg.eval.costQ;
evalOptions.needDetail = false;
evalOptions.enablePruning = cfg.eval.enablePruning;
evalOptions.incumbentCost = incumbentCost;
evalOptions.method = cfg.eval.method;

evalResult = case1_evaluate_scheme(dataset, scheme, evalOptions);
costValue = evalResult.meanCost;
processedSamples = evalResult.processedSampleCount;
wasPruned = evalResult.pruned;
exactFlag = false;

% Search scoring deliberately uses a representative subset. Full-train
% records are created only by exact_evaluate_scheme_SUBMIT for generation elites.
end

function [generationBestExact, globalBest, bestByCount, archive, exactCache] = exact_update_candidates_SUBMIT( ...
        candidateMatrix, dataset, sampleOrder, cfg, priority, globalBest, bestByCount, archive, exactCache)
generationBestExact = globalBest.cost;
exactRecords = repmat(empty_record_SUBMIT(), 0, 1);

for rowId = 1:size(candidateMatrix, 1)
    [record, exactCache] = exact_evaluate_scheme_SUBMIT(candidateMatrix(rowId, :), dataset, sampleOrder, cfg, priority, exactCache);
    exactRecords(end + 1, 1) = record; %#ok<AGROW>
end

if ~isempty(exactRecords)
    generationBestExact = min([exactRecords.cost]);
    [globalBest, bestByCount, archive] = assimilate_exact_records_SUBMIT(exactRecords, globalBest, bestByCount, archive, cfg);
end
end

function [record, exactCache] = exact_evaluate_scheme_SUBMIT(scheme, dataset, sampleOrder, cfg, priority, exactCache)
scheme = repair_binary_scheme_SUBMIT(scheme, cfg.search, priority);
key = scheme_key_SUBMIT(scheme);

if isKey(exactCache, key)
    record = exactCache(key);
    return;
end

evalOptions.sampleIndices = 1:dataset.sampleCount;
evalOptions.sampleOrder = sampleOrder;
evalOptions.costQ = cfg.eval.costQ;
evalOptions.needDetail = false;
evalOptions.enablePruning = false;
evalOptions.incumbentCost = inf;
evalOptions.method = cfg.eval.method;

evalResult = case1_evaluate_scheme(dataset, scheme, evalOptions);
record = record_from_eval_result_SUBMIT(scheme, evalResult);
exactCache(key) = record;
end

function [globalBest, bestByCount, archive] = assimilate_exact_records_SUBMIT(records, globalBest, bestByCount, archive, cfg)
if isempty(records)
    return;
end

if ~isstruct(records)
    return;
end

[~, order] = sort([records.cost], 'ascend');
records = records(order);

for i = 1:numel(records)
    record = records(i);
    if isempty(record.scheme)
        continue;
    end

    slot = record.count - cfg.search.minPointCount + 1;
    if slot >= 1 && slot <= numel(bestByCount) && record.cost < bestByCount(slot).cost
        bestByCount(slot) = record;
    end

    if record.cost < globalBest.cost
        globalBest = record;
    end
end

archive = archive_update_SUBMIT(archive, records, cfg.search.archiveSize, cfg.search.archiveMinDistance);
end

function archive = archive_update_SUBMIT(archive, records, archiveSize, minDistance)
if isempty(records)
    return;
end

for i = 1:numel(records)
    record = records(i);
    if isempty(record.scheme)
        continue;
    end

    replaceId = 0;
    bestDistance = inf;
    for arcId = 1:numel(archive)
        distanceValue = archive_distance_SUBMIT(record.scheme, archive(arcId).scheme);
        if distanceValue < bestDistance
            bestDistance = distanceValue;
            replaceId = arcId;
        end
    end

    if isempty(archive)
        archive = record;
    elseif bestDistance < minDistance
        if record.cost < archive(replaceId).cost
            archive(replaceId) = record;
        end
    else
        archive(end + 1) = record; %#ok<AGROW>
    end
end

[~, order] = sort([archive.cost], 'ascend');
archive = archive(order);
if numel(archive) > archiveSize
    archive = archive(1:archiveSize);
end
end

function distanceValue = archive_distance_SUBMIT(scheme1, scheme2)
distanceValue = nnz(xor(scheme1 > 0.5, scheme2 > 0.5));
end

function adaptive = adaptive_params_SUBMIT(search, gen, maxGenerations, stagnationCount, islands)
progress = (gen - 1) / max(1, maxGenerations - 1);
adaptive.mutationProb = search.initialMutationProb * (1 - progress) + search.finalMutationProb * progress;
adaptive.destroyRepairProb = search.destroyRepairProb;
adaptive.moveRadius = 1 + floor(3 * progress);
adaptive.countJitter = search.countJitter;
adaptive.immigrantCount = search.immigrantCount;

diversityList = zeros(numel(islands), 1);
for islandId = 1:numel(islands)
    diversityList(islandId) = compute_population_diversity_SUBMIT(islands{islandId});
end

if mean(diversityList) < search.lowDiversityThreshold
    adaptive.mutationProb = min(0.55, adaptive.mutationProb + 0.08);
    adaptive.destroyRepairProb = min(0.50, adaptive.destroyRepairProb + 0.08);
end

if stagnationCount >= search.stagnationWindow
    adaptive.mutationProb = min(0.60, adaptive.mutationProb + search.stagnationMutationBoost);
    adaptive.destroyRepairProb = min(0.55, adaptive.destroyRepairProb + search.stagnationDestroyBoost);
    adaptive.immigrantCount = max(search.immigrantCount + 2, ceil(1.5 * search.immigrantCount));
end
end

function islands = reproduce_islands_SUBMIT(islandSummaries, dataset, cfg, priority, archive, adaptive)
search = cfg.search;
geneCount = dataset.pointCount;
popSize = search.populationSizePerIsland;
islands = cell(search.islandCount, 1);

for islandId = 1:search.islandCount
    summary = islandSummaries(islandId);
    population = summary.population;
    searchCost = summary.searchCost;

    nextPopulation = zeros(popSize, geneCount);
    eliteCount = min(search.eliteCount, popSize);
    nextPopulation(1:eliteCount, :) = population(1:eliteCount, :);

    immigrantCount = min(adaptive.immigrantCount, max(0, popSize - eliteCount));
    if immigrantCount > 0
        nextPopulation(eliteCount + 1:eliteCount + immigrantCount, :) = ...
            random_binary_population_SUBMIT(immigrantCount, geneCount, search, priority);
    end

    cursor = eliteCount + immigrantCount + 1;
    while cursor <= popSize
        parent1 = population(tournament_pick_SUBMIT(searchCost, search.tournamentSize), :);
        parent2 = population(tournament_pick_SUBMIT(searchCost, search.tournamentSize), :);

        if rand < search.crossoverProb
            child1 = crossover_binary_SUBMIT(parent1, parent2, search, priority);
            child2 = crossover_binary_SUBMIT(parent2, parent1, search, priority);
        else
            child1 = parent1;
            child2 = parent2;
        end

        child1 = mutate_binary_scheme_SUBMIT(child1, geneCount, search, priority, adaptive.mutationProb, adaptive.moveRadius);
        child2 = mutate_binary_scheme_SUBMIT(child2, geneCount, search, priority, adaptive.mutationProb, adaptive.moveRadius);

        if rand < adaptive.destroyRepairProb
            child1 = destroy_repair_SUBMIT(child1, geneCount, search, priority);
        end
        if rand < adaptive.destroyRepairProb
            child2 = destroy_repair_SUBMIT(child2, geneCount, search, priority);
        end

        nextPopulation(cursor, :) = child1;
        if cursor + 1 <= popSize
            nextPopulation(cursor + 1, :) = child2;
        end
        cursor = cursor + 2;
    end

    nextPopulation = inject_archive_children_SUBMIT(nextPopulation, archive, search, priority);
    nextPopulation = enforce_population_diversity_SUBMIT(nextPopulation, search, priority);
    islands{islandId} = nextPopulation;
end
end

function child = crossover_binary_SUBMIT(parent1, parent2, search, priority)
geneCount = numel(parent1);
count1 = nnz(parent1);
count2 = nnz(parent2);
targetCount = round((count1 + count2) / 2 + randi([-search.countJitter, search.countJitter], 1, 1));
targetCount = min(search.maxPointCount, max(search.minPointCount, targetCount));

modeId = randi(3);
switch modeId
    case 1
        mask = rand(1, geneCount) < 0.5;
        child = parent1;
        child(mask) = parent2(mask);

    case 2
        child = double(parent1 & parent2);
        unionIdx = find((parent1 | parent2) & ~child);
        if ~isempty(unionIdx)
            score = priority.score(unionIdx) + 0.25 * (parent1(unionIdx) + parent2(unionIdx));
            pickCount = min(numel(unionIdx), max(0, targetCount - nnz(child)));
            chosen = weighted_sample_without_replacement_SUBMIT(unionIdx, score, pickCount);
            child(chosen) = 1;
        end

    otherwise
        cut = sort(randperm(geneCount - 1, 2));
        child = [parent1(1:cut(1)), parent2(cut(1) + 1:cut(2)), parent1(cut(2) + 1:end)];
end

child = repair_binary_scheme_SUBMIT(child, search, priority);
end

function child = mutate_binary_scheme_SUBMIT(child, geneCount, search, priority, mutationProb, moveRadius)
child = double(child > 0.5);

if rand < mutationProb
    selectedIdx = find(child);
    if ~isempty(selectedIdx)
        oldPos = selectedIdx(randi(numel(selectedIdx)));
        newPos = min(geneCount, max(1, oldPos + randi([-moveRadius, moveRadius], 1, 1)));
        child(oldPos) = 0;
        child(newPos) = 1;
    end
end

if rand < search.addProb * mutationProb
    pool = build_add_pool_SUBMIT(child, priority, search, 8);
    if ~isempty(pool)
        child(pool(randi(numel(pool)))) = 1;
    end
end

if rand < search.deleteProb * mutationProb
    selectedIdx = find(child);
    if numel(selectedIdx) > search.minPointCount
        removeScore = compute_removal_score_SUBMIT(selectedIdx, priority.score);
        [~, order] = sort(removeScore, 'descend');
        child(selectedIdx(order(1))) = 0;
    end
end

if rand < search.swapProb * mutationProb
    selectedIdx = find(child);
    pool = build_add_pool_SUBMIT(child, priority, search, 6);
    if ~isempty(selectedIdx) && ~isempty(pool)
        removeScore = compute_removal_score_SUBMIT(selectedIdx, priority.score);
        [~, order] = sort(removeScore, 'descend');
        child(selectedIdx(order(1))) = 0;
        child(pool(randi(numel(pool)))) = 1;
    end
end

child = repair_binary_scheme_SUBMIT(child, search, priority);
end

function child = destroy_repair_SUBMIT(child, geneCount, search, priority)
child = double(child > 0.5);
selectedIdx = find(child);
if isempty(selectedIdx)
    child = random_binary_scheme_SUBMIT(geneCount, search, priority);
    return;
end

removeCount = randi([search.destroyRepairRemoveMin, search.destroyRepairRemoveMax], 1, 1);
removeCount = min(removeCount, max(1, numel(selectedIdx) - search.minPointCount + 1));

removeScore = compute_removal_score_SUBMIT(selectedIdx, priority.score);
[~, order] = sort(removeScore, 'descend');
toRemove = selectedIdx(order(1:removeCount));
child(toRemove) = 0;

pool = build_add_pool_SUBMIT(child, priority, search, 12);
addCount = min(removeCount, numel(pool));
if addCount > 0
    chosen = weighted_sample_without_replacement_SUBMIT(pool, priority.score(pool), addCount);
    child(chosen) = 1;
end

if nnz(child) < search.minPointCount
    child = repair_binary_scheme_SUBMIT(child, search, priority);
end
end

function pool = build_add_pool_SUBMIT(scheme, priority, search, maxCount)
selectedIdx = find(scheme);
geneCount = numel(scheme);
candidate = zeros(1, 0);

if ~isempty(selectedIdx)
    edges = [0, selectedIdx, geneCount + 1];
    gapLength = diff(edges) - 1;
    [~, gapOrder] = sort(gapLength, 'descend');
    for id = 1:numel(gapOrder)
        gapId = gapOrder(id);
        left = edges(gapId) + 1;
        right = edges(gapId + 1) - 1;
        if left <= right
            mid = round((left + right) / 2);
            candidate(end + 1) = mid; %#ok<AGROW>
            if left < mid
                candidate(end + 1) = left; %#ok<AGROW>
            end
            if right > mid
                candidate(end + 1) = right; %#ok<AGROW>
            end
        end
        if numel(candidate) >= 3 * maxCount
            break;
        end
    end

    for pos = selectedIdx
        if pos > 1
            candidate(end + 1) = pos - 1; %#ok<AGROW>
        end
        if pos < geneCount
            candidate(end + 1) = pos + 1; %#ok<AGROW>
        end
    end
end

candidate = [candidate, priority.order(1:min(numel(priority.order), 3 * maxCount))]; %#ok<AGROW>
candidate = unique(candidate(candidate >= 1 & candidate <= geneCount), 'stable');
candidate = candidate(~scheme(candidate));

if isempty(candidate)
    pool = [];
    return;
end

[~, order] = sort(priority.score(candidate), 'descend');
pool = candidate(order(1:min(maxCount, numel(order))));
end

function population = inject_archive_children_SUBMIT(population, archive, search, priority)
if isempty(archive) || mod(size(population, 1), max(1, search.archiveInjectEvery)) ~= 0
    % Keep the condition simple and deterministic to avoid excessive churn.
end

injectCount = min(search.archiveInjectCount, numel(archive));
for injectId = 1:injectCount
    parent1 = archive(injectId).scheme;
    parent2 = archive(randi(numel(archive))).scheme;
    child = crossover_binary_SUBMIT(parent1, parent2, search, priority);
    child = mutate_binary_scheme_SUBMIT(child, numel(child), search, priority, 0.20, 3);
    population(end - injectId + 1, :) = child;
end
end

function islands = migrate_islands_SUBMIT(islands, islandSummaries, migrationCount)
if migrationCount <= 0
    return;
end

islandCount = numel(islands);
for islandId = 1:islandCount
    source = islandSummaries(islandId).population(1:min(migrationCount, size(islands{islandId}, 1)), :);
    targetId = mod(islandId, islandCount) + 1;
    replaceCount = size(source, 1);
    islands{targetId}(end - replaceCount + 1:end, :) = source;
end
end

function population = enforce_population_diversity_SUBMIT(population, search, priority)
population = unique(population, 'rows', 'stable');
targetSize = search.populationSizePerIsland;
geneCount = size(population, 2);

while size(population, 1) < targetSize
    newcomer = random_binary_population_SUBMIT(1, geneCount, search, priority);
    population = unique([population; newcomer], 'rows', 'stable'); %#ok<AGROW>
end

population = population(1:targetSize, :);
end

function population = random_binary_population_SUBMIT(count, geneCount, search, priority)
population = zeros(count, geneCount);
for rowId = 1:count
    population(rowId, :) = random_binary_scheme_SUBMIT(geneCount, search, priority);
end
end

function scheme = random_binary_scheme_SUBMIT(geneCount, search, priority)
rawCount = round(search.countSoftCenter + randn() * 2.2);
targetCount = min(search.maxPointCount, max(search.minPointCount, rawCount));
candidate = weighted_sample_without_replacement_SUBMIT(1:geneCount, priority.score, targetCount);
scheme = scheme_from_indices_SUBMIT(candidate, geneCount);
scheme = repair_binary_scheme_SUBMIT(scheme, search, priority);
end

function [bestRecord, exactCache] = refine_scheme_variable_SUBMIT(dataset, baseRecord, sampleOrder, cfg, priority, exactCache, options)
bestRecord = baseRecord;
if isempty(bestRecord.scheme)
    return;
end

for passId = 1:options.maxPasses
    improved = false;
    currentScheme = bestRecord.scheme;
    selectedIdx = find(currentScheme);
    addPool = build_add_pool_SUBMIT(currentScheme, priority, cfg.search, options.addCandidateCount);

    if numel(selectedIdx) > cfg.search.minPointCount
        removeScore = compute_removal_score_SUBMIT(selectedIdx, priority.score);
        [~, removeOrder] = sort(removeScore, 'descend');
        removeCandidates = selectedIdx(removeOrder(1:min(options.removeCandidateCount, numel(removeOrder))));
        for removeId = 1:numel(removeCandidates)
            candidate = currentScheme;
            candidate(removeCandidates(removeId)) = 0;
            [bestRecord, exactCache, improved] = try_candidate_SUBMIT( ...
                candidate, dataset, sampleOrder, cfg, priority, exactCache, bestRecord, improved);
        end
    end

    if numel(selectedIdx) < cfg.search.maxPointCount
        for addId = 1:numel(addPool)
            candidate = currentScheme;
            candidate(addPool(addId)) = 1;
            [bestRecord, exactCache, improved] = try_candidate_SUBMIT( ...
                candidate, dataset, sampleOrder, cfg, priority, exactCache, bestRecord, improved);
        end
    end

    removeScore = compute_removal_score_SUBMIT(selectedIdx, priority.score);
    [~, removeOrder] = sort(removeScore, 'descend');
    removeCandidates = selectedIdx(removeOrder(1:min(options.swapRemoveCandidateCount, numel(removeOrder))));
    swapAddPool = build_add_pool_SUBMIT(currentScheme, priority, cfg.search, options.swapAddCandidateCount);

    for removeId = 1:numel(removeCandidates)
        baseCandidate = currentScheme;
        baseCandidate(removeCandidates(removeId)) = 0;

        shiftRange = max(1, removeCandidates(removeId) - options.moveWindow): ...
            min(numel(currentScheme), removeCandidates(removeId) + options.moveWindow);
        for newPos = shiftRange
            moveCandidate = baseCandidate;
            moveCandidate(newPos) = 1;
            [bestRecord, exactCache, improved] = try_candidate_SUBMIT( ...
                moveCandidate, dataset, sampleOrder, cfg, priority, exactCache, bestRecord, improved);
        end

        for addId = 1:numel(swapAddPool)
            swapCandidate = baseCandidate;
            swapCandidate(swapAddPool(addId)) = 1;
            [bestRecord, exactCache, improved] = try_candidate_SUBMIT( ...
                swapCandidate, dataset, sampleOrder, cfg, priority, exactCache, bestRecord, improved);
        end
    end

    if ~improved
        break;
    end
end
end

function [bestRecord, exactCache, improved] = try_candidate_SUBMIT( ...
        candidate, dataset, sampleOrder, cfg, priority, exactCache, bestRecord, improved)
[record, exactCache] = exact_evaluate_scheme_SUBMIT(candidate, dataset, sampleOrder, cfg, priority, exactCache);
if record.cost + 1e-9 < bestRecord.cost
    bestRecord = record;
    improved = true;
end
end

function record = record_from_eval_result_SUBMIT(scheme, evalResult)
record = empty_record_SUBMIT();
record.scheme = scheme(:)';
record.cost = evalResult.meanCost;
record.count = evalResult.selectedCount;
record.selectedIdx = evalResult.selectedIdx(:)';
record.temperature = evalResult.selectedTemperature(:)';
record.sampleCost = evalResult.sampleCost(:);
end

function record = empty_record_SUBMIT()
record.scheme = [];
record.cost = inf;
record.count = NaN;
record.selectedIdx = [];
record.temperature = [];
record.sampleCost = [];
end

function score = compute_removal_score_SUBMIT(selectedIdx, priorityScore)
if isempty(selectedIdx)
    score = [];
    return;
end

leftGap = [inf, diff(selectedIdx)];
rightGap = [diff(selectedIdx), inf];
localGap = min(leftGap, rightGap);
importance = priorityScore(selectedIdx);
importance = normalize01_SUBMIT(importance);
compactness = 1 ./ (1 + localGap);
score = compactness + (1 - importance);
end

function order = build_priority_order_SUBMIT(fullSampleIndices, hardSamplePool)
if isempty(hardSamplePool)
    order = fullSampleIndices;
    return;
end

hardSamplePool = hardSamplePool(:)';
hardSamplePool = hardSamplePool(ismember(hardSamplePool, fullSampleIndices));
missing = setdiff(fullSampleIndices, hardSamplePool, 'stable');
order = [hardSamplePool, missing];
end

function sampleIndices = build_search_sample_indices_SUBMIT(sampleCount, targetCount, randomCount, seed)
targetCount = min(targetCount, sampleCount);
gridCount = max(0, targetCount - randomCount);
gridIndices = unique(round(linspace(1, sampleCount, max(1, gridCount))));

rng(seed + 7919, 'twister');
randomCount = min(randomCount, sampleCount);
randomIndices = randperm(sampleCount, randomCount);
sampleIndices = unique([gridIndices, randomIndices], 'stable');

if numel(sampleIndices) > targetCount
    sampleIndices = sampleIndices(1:targetCount);
end
end

function sampleIndices = build_iteration_sample_indices_SUBMIT(baseIndices, hardSamplePool, fullSampleIndices, targetCount)
if isempty(hardSamplePool)
    sampleIndices = baseIndices;
else
    hardSamplePool = hardSamplePool(:)';
    hardSamplePool = hardSamplePool(ismember(hardSamplePool, fullSampleIndices));
    sampleIndices = unique([hardSamplePool, baseIndices], 'stable');
end

if numel(sampleIndices) > targetCount
    keepHard = hardSamplePool(ismember(hardSamplePool, sampleIndices));
    remaining = setdiff(sampleIndices, keepHard, 'stable');
    remaining = remaining(1:max(0, targetCount - numel(keepHard)));
    sampleIndices = unique([keepHard, remaining], 'stable');
end
end

function diversity = compute_population_diversity_SUBMIT(population)
if size(population, 1) <= 1
    diversity = 0;
    return;
end

sampleSize = min(12, size(population, 1));
pick = randperm(size(population, 1), sampleSize);
pairCount = 0;
distanceSum = 0;

for i = 1:sampleSize
    for j = i + 1:sampleSize
        pairCount = pairCount + 1;
        distanceSum = distanceSum + nnz(xor(population(pick(i), :), population(pick(j), :)));
    end
end

diversity = distanceSum / max(1, pairCount) / size(population, 2);
end

function idx = weighted_sample_without_replacement_SUBMIT(pool, weights, pickCount)
pool = pool(:)';
weights = weights(:)';
pickCount = min(pickCount, numel(pool));

if pickCount <= 0 || isempty(pool)
    idx = [];
    return;
end

weights = max(weights, eps);
idx = zeros(1, pickCount);
localPool = pool;
localWeights = weights;

for pickId = 1:pickCount
    cdf = cumsum(localWeights) / sum(localWeights);
    pos = find(cdf >= rand, 1, 'first');
    idx(pickId) = localPool(pos);
    localPool(pos) = [];
    localWeights(pos) = [];
    if isempty(localPool)
        idx = idx(1:pickId);
        break;
    end
end
end

function scheme = repair_binary_scheme_SUBMIT(rawScheme, search, priority)
scheme = double(rawScheme(:)' > 0.5);
geneCount = numel(scheme);
selectedCount = nnz(scheme);

if selectedCount < search.minPointCount
    need = search.minPointCount - selectedCount;
    complement = find(~scheme);
    chosen = weighted_sample_without_replacement_SUBMIT(complement, priority.score(complement), need);
    scheme(chosen) = 1;
elseif selectedCount > search.maxPointCount
    selectedIdx = find(scheme);
    removeScore = compute_removal_score_SUBMIT(selectedIdx, priority.score);
    [~, order] = sort(removeScore, 'descend');
    dropCount = selectedCount - search.maxPointCount;
    scheme(selectedIdx(order(1:dropCount))) = 0;
end
end

function idx = indices_from_scheme_SUBMIT(scheme)
idx = find(scheme > 0.5);
end

function scheme = scheme_from_indices_SUBMIT(indices, geneCount)
scheme = zeros(1, geneCount);
indices = unique(round(indices(:)'));
indices = indices(indices >= 1 & indices <= geneCount);
scheme(indices) = 1;
end

function pos = tournament_pick_SUBMIT(costVector, tournamentSize)
pick = randperm(numel(costVector), min(tournamentSize, numel(costVector)));
[~, localId] = min(costVector(pick));
pos = pick(localId);
end

function key = scheme_key_SUBMIT(scheme)
key = char(scheme(:)' + '0');
end

function y = normalize01_SUBMIT(x)
x = x(:)';
xMin = min(x);
xMax = max(x);
if xMax <= xMin + eps
    y = zeros(size(x));
else
    y = (x - xMin) / (xMax - xMin);
end
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


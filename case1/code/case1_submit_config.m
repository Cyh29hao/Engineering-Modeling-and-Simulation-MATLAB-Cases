function cfg = case1_submit_config(profile)
%case1_submit_config Configuration for variable-point-count Case1 search.

if nargin < 1 || isempty(profile)
    profile = 'balanced';
end

cfg.meta.profile = profile;
cfg.meta.versionName = 'Case1 submit';

cfg.search.minPointCount = 4;
cfg.search.maxPointCount = 14;
cfg.search.populationSizePerIsland = 30;
cfg.search.islandCount = 4;
cfg.search.maxGenerations = 42;
cfg.search.eliteCount = 2;
cfg.search.immigrantCount = 3;
cfg.search.exactEliteCount = 4;
cfg.search.tournamentSize = 4;
cfg.search.crossoverProb = 0.90;
cfg.search.initialMutationProb = 0.24;
cfg.search.finalMutationProb = 0.12;
cfg.search.addProb = 0.18;
cfg.search.deleteProb = 0.16;
cfg.search.moveProb = 0.34;
cfg.search.swapProb = 0.22;
cfg.search.destroyRepairProb = 0.18;
cfg.search.destroyRepairRemoveMin = 1;
cfg.search.destroyRepairRemoveMax = 3;
cfg.search.seedJitterCount = 18;
cfg.search.seedJitterRadius = 6;
cfg.search.archiveSize = 40;
cfg.search.archiveMinDistance = 8;
cfg.search.archiveInjectEvery = 4;
cfg.search.archiveInjectCount = 2;
cfg.search.migrationEvery = 6;
cfg.search.migrationCount = 2;
cfg.search.lowDiversityThreshold = 0.22;
cfg.search.stagnationWindow = 8;
cfg.search.stagnationMutationBoost = 0.10;
cfg.search.stagnationDestroyBoost = 0.12;
cfg.search.countJitter = 2;
cfg.search.countSoftCenter = 7;
cfg.search.countSoftPenalty = 0.02;

cfg.eval.costQ = 70;
cfg.eval.method = 'spline';
cfg.eval.enablePruning = true;
cfg.eval.searchSampleCount = 360;
cfg.eval.searchRandomSampleCount = 120;
cfg.eval.hardSampleCount = 48;

cfg.localSearch.enable = true;
cfg.localSearch.every = 4;
cfg.localSearch.maxPasses = 1;
cfg.localSearch.moveWindow = 2;
cfg.localSearch.addCandidateCount = 12;
cfg.localSearch.removeCandidateCount = 6;
cfg.localSearch.swapAddCandidateCount = 6;
cfg.localSearch.swapRemoveCandidateCount = 4;
cfg.localSearch.finalPasses = 3;
cfg.localSearch.finalMoveWindow = 4;
cfg.localSearch.finalAddCandidateCount = 18;
cfg.localSearch.finalRemoveCandidateCount = 8;

cfg.random.baseSeed = 20260423;
cfg.random.restartCount = 2;

cfg.finalRefine.enableTwoDeleteTwoAdd = true;
cfg.finalRefine.enableThreeDeleteThreeAdd = false;
cfg.finalRefine.enableCountTransition = true;
cfg.finalRefine.countTransitionTopCount = 3;
cfg.finalRefine.countTransitionMaxSourceCount = 8;
cfg.finalRefine.countTransitionWindow = 10;
cfg.finalRefine.countTransitionPasses = 4;
cfg.finalRefine.twoDeleteSourceCount = 3;
cfg.finalRefine.twoDeleteMaxStartCost = 500;
cfg.finalRefine.twoDeletePasses = 3;

cfg.output.consoleEvery = 1;
cfg.output.checkpointEvery = 0;
cfg.output.writeHistoryCsv = true;
cfg.output.writeSummaryTxt = true;

switch lower(profile)
    case 'balanced'
        % Keep defaults.

    case 'aggressive'
        cfg.search.populationSizePerIsland = 36;
        cfg.search.maxGenerations = 56;
        cfg.search.exactEliteCount = 5;
        cfg.search.seedJitterCount = 24;
        cfg.search.archiveSize = 56;
        cfg.search.archiveInjectCount = 3;
        cfg.search.stagnationMutationBoost = 0.14;
        cfg.search.stagnationDestroyBoost = 0.16;
        cfg.localSearch.finalPasses = 4;
        cfg.localSearch.finalMoveWindow = 5;
        cfg.random.restartCount = 3;

    case 'fast'
        cfg.search.maxPointCount = 12;
        cfg.search.populationSizePerIsland = 20;
        cfg.search.islandCount = 3;
        cfg.search.maxGenerations = 18;
        cfg.search.exactEliteCount = 3;
        cfg.search.seedJitterCount = 4;
        cfg.search.archiveSize = 10;
        cfg.search.archiveInjectEvery = 2;
        cfg.search.archiveInjectCount = 1;
        cfg.search.migrationEvery = 2;
        cfg.search.migrationCount = 1;
        cfg.search.stagnationWindow = 3;
        cfg.eval.searchSampleCount = 260;
        cfg.eval.searchRandomSampleCount = 80;
        cfg.eval.hardSampleCount = 32;
        cfg.localSearch.every = 3;
        cfg.localSearch.finalPasses = 2;
        cfg.random.restartCount = 1;
        cfg.output.checkpointEvery = 3;

    case 'submit'
        cfg.search.maxPointCount = 12;
        cfg.search.populationSizePerIsland = 20;
        cfg.search.islandCount = 1;
        cfg.search.maxGenerations = 5;
        cfg.search.eliteCount = 2;
        cfg.search.immigrantCount = 2;
        cfg.search.exactEliteCount = 2;
        cfg.search.seedJitterCount = 4;
        cfg.search.archiveSize = 8;
        cfg.search.archiveInjectEvery = 2;
        cfg.search.archiveInjectCount = 1;
        cfg.search.migrationEvery = 100;
        cfg.search.migrationCount = 0;
        cfg.search.stagnationWindow = 3;
        cfg.eval.searchSampleCount = 120;
        cfg.eval.searchRandomSampleCount = 40;
        cfg.eval.hardSampleCount = 16;
        cfg.localSearch.enable = true;
        cfg.localSearch.every = 1;
        cfg.localSearch.maxPasses = 1;
        cfg.localSearch.moveWindow = 3;
        cfg.localSearch.addCandidateCount = 8;
        cfg.localSearch.removeCandidateCount = 5;
        cfg.localSearch.swapAddCandidateCount = 5;
        cfg.localSearch.swapRemoveCandidateCount = 4;
        cfg.localSearch.finalPasses = 1;
        cfg.localSearch.finalMoveWindow = 3;
        cfg.localSearch.finalAddCandidateCount = 18;
        cfg.localSearch.finalRemoveCandidateCount = 8;
        cfg.random.restartCount = 1;
        cfg.finalRefine.enableTwoDeleteTwoAdd = false;
        cfg.finalRefine.enableThreeDeleteThreeAdd = false;
        cfg.finalRefine.enableCountTransition = false;
        cfg.output.checkpointEvery = 0;

    case 'smoke'
        cfg.search.maxPointCount = 10;
        cfg.search.populationSizePerIsland = 8;
        cfg.search.islandCount = 2;
        cfg.search.maxGenerations = 4;
        cfg.search.eliteCount = 1;
        cfg.search.immigrantCount = 1;
        cfg.search.exactEliteCount = 2;
        cfg.search.seedJitterCount = 4;
        cfg.search.archiveSize = 10;
        cfg.search.archiveInjectEvery = 2;
        cfg.search.archiveInjectCount = 1;
        cfg.search.migrationEvery = 2;
        cfg.search.migrationCount = 1;
        cfg.search.stagnationWindow = 3;
        cfg.eval.hardSampleCount = 8;
        cfg.eval.searchSampleCount = 80;
        cfg.eval.searchRandomSampleCount = 24;
        cfg.localSearch.every = 2;
        cfg.localSearch.maxPasses = 1;
        cfg.localSearch.moveWindow = 1;
        cfg.localSearch.addCandidateCount = 6;
        cfg.localSearch.removeCandidateCount = 4;
        cfg.localSearch.swapAddCandidateCount = 4;
        cfg.localSearch.swapRemoveCandidateCount = 3;
        cfg.localSearch.finalPasses = 1;
        cfg.localSearch.finalMoveWindow = 2;
        cfg.localSearch.finalAddCandidateCount = 8;
        cfg.localSearch.finalRemoveCandidateCount = 4;
        cfg.random.restartCount = 1;
        cfg.finalRefine.enableTwoDeleteTwoAdd = false;
        cfg.finalRefine.enableThreeDeleteThreeAdd = false;
        cfg.finalRefine.enableCountTransition = false;

    otherwise
        error('Unknown profile: %s', profile);
end
end


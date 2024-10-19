function [keyA, keyB, bestScore, xcorrAtBest, flipTFAtBest, circShiftAtBest, cache, cacheKeys] = get_best_pair_indices(consensusKeyPool, barcodeStructsMap, cache)
    barcodeKeys = consensusKeyPool;
    ignoreScoresMat = tril(true(length(barcodeKeys))); % nonsensical pairings against self and redundancies (where the pairing order/reorientation strategy is just in reverse but score-wise results should be duplicates due to symmetrical scoring)

    import CBT.Consensus.CoreCacheBridge.get_pairwise_similarities;
    [pairwiseBsos, cache] = get_pairwise_similarities(barcodeKeys, barcodeStructsMap, cache, true);
    pairwiseBestScores = cellfun(@(bsosStruct) bsosStruct.bestScore, pairwiseBsos);

    pairwiseBestScores(ignoreScoresMat) = NaN;
    bestPairwiseBestScore = max(pairwiseBestScores(:)); % note that this ignores NaNs
    [keyIdxA, keyIdxB] = ind2sub(size(pairwiseBestScores), find(pairwiseBestScores == bestPairwiseBestScore, 1, 'first'));

    bestBsosStruct = pairwiseBsos{keyIdxA, keyIdxB};
    keyA = barcodeKeys{keyIdxA};
    keyB = barcodeKeys{keyIdxB};

    bestScore = bestBsosStruct.bestScore;
    xcorrAtBest = bestBsosStruct.xcorrAtBest;
    flipTFAtBest = bestBsosStruct.flipTFAtBest;
    circShiftAtBest = bestBsosStruct.circShiftAtBest;

    cacheKeys = cellfun(@(bsosStruct) bsosStruct.cacheKey, pairwiseBsos, 'UniformOutput', false);
    cacheKeys = cacheKeys(:);
    cacheKeys = cacheKeys(~cellfun(@isempty, cacheKeys));
end
function [isInCache, cacheKey, cachedValue] = check_cache_for_best_synced_orientation_similarity(circSeqA, circSeqB, indexWeightsA, indexWeightsB, cache)
    import Fancy.Utils.data_hash;

    defaultUncachedValue = [];
    cachedValue = defaultUncachedValue;

    cacheKey = data_hash({circSeqA, circSeqB, indexWeightsA, indexWeightsB});

    isInCache = isKey(cache, cacheKey);
    if isInCache
        cachedValue = cache(cacheKey);
    end

    % % Todo(?) note: if we wanted to be clever, for circular
    % % sequences of the same length where the reverse ordering of
    % % the pair of sequences+weights has results for sync orientating
    % % cached, we should theoretically be able to quickly compute 
    % % some valid results using something like the code below.
    % % In the unlikely case that there is a tie for the best score
    % % however, the orientation results could theoretically differ
    % % than if computed directly.
    % % (Currently choosing not to do this because it also
    % % might make the code a bit more complicated/fragile.)
    %
    % if not(isInCache)
    %     reverseOrderCacheKey = data_hash({circSeqB, circSeqA, indexWeightsB, indexWeightsA});
    %     isInCache = isKey(cache, reverseOrderCacheKey);
    %     reverseOrderCachedValue = cache(reverseOrderCacheKey);
    %     seqLen = length(circSeqB);
    %     if not(reverseOrderCachedValue.flipTFAtBest)
    %         reverseOrderCachedValue.circShiftAtBest = (-1 * reverseOrderCachedValue.circShiftAtBest) - seqLen;
    %     end
    % end
end
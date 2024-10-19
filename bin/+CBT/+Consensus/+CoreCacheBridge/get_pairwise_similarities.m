function [pairwiseBsos, cache] = get_pairwise_similarities(barcodeKeys, barcodeStructsMap, cache, allowRecalculate)
    validateattributes(barcodeKeys, {'cell'}, {'column'}, 1)

    validateattributes(barcodeStructsMap, {'containers.Map'}, {}, 2);
    if not(all(cellfun(@(barcodeKey) isKey(barcodeStructsMap, barcodeKey), barcodeKeys)))
        error('Missing barcodes for some barcode keys');
    end

    if (nargin < 3) || isempty(cache)
        cache = containers.Map();
    else
        validateattributes(cache, {'containers.Map'}, {}, 3);
    end
    if (nargin < 4)
        allowRecalculate = false;
    else
        validateattributes(allowRecalculate, {'logical'}, {'scalar'}, 4);
    end

    barcodeLens = cellfun(@(barcodeKey) length(barcodeStructsMap(barcodeKey).barcode), barcodeKeys);
    barcodeIndexWeightLens = cellfun(@(barcodeKey) length(barcodeStructsMap(barcodeKey).indexWeights), barcodeKeys);

    if any(diff([barcodeLens; barcodeIndexWeightLens]) ~= 0)
        error('Length inconsistency with barcodes and/or weights');
    end

    numBarcodes = numel(barcodeKeys);
    pairwiseBsos = cell(numBarcodes, numBarcodes);

    import CBT.Consensus.Caching.check_cache_for_best_synced_orientation_similarity;
    import CBT.Consensus.Core.calc_best_synced_orientation_similarity;

    for barcodeIdxA = 1:numBarcodes

        barcodeKeyA = barcodeKeys{barcodeIdxA};
        barcodeStructA = barcodeStructsMap(barcodeKeyA);
        barcodeA = barcodeStructA.barcode;
        indexWeightsA = barcodeStructA.indexWeights;

        for barcodeIdxB = 1:numBarcodes

            barcodeKeyB = barcodeKeys{barcodeIdxB};
            barcodeStructB = barcodeStructsMap(barcodeKeyB);
            barcodeB = barcodeStructB.barcode;
            indexWeightsB = barcodeStructB.indexWeights;


            [isInCache, cacheKey, bsosStruct] = check_cache_for_best_synced_orientation_similarity(barcodeA, barcodeB, indexWeightsA, indexWeightsB, cache);
            if not(isInCache)
                if allowRecalculate
                    bsosStruct = calc_best_synced_orientation_similarity(barcodeA, barcodeB, indexWeightsA, indexWeightsB);
                    cache(cacheKey) = bsosStruct;
                else
                    bsosStruct = struct();
                end
            end
            bsosStruct.cacheKey = cacheKey;

            pairwiseBsos{barcodeIdxA, barcodeIdxB} = bsosStruct;

        end
    end
end
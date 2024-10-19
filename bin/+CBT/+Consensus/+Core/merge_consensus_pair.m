function [circSeqAB, indexWeightsAB] = merge_consensus_pair(circSeqA, circSeqB, indexWeightsA, indexWeightsB, flipTFAtBest, circShiftAtBest)
    import Barcoding.Reorienting.reorient_barcode_and_bitmask;

    if length(circSeqA) ~= length(circSeqB)
        error('Sequences should be the same length!');
    end

    [syncedCircSeqB, syncedIndexWeightsB] = reorient_barcode_and_bitmask(circSeqB, indexWeightsB, flipTFAtBest, circShiftAtBest);

    circSeqAB = circSeqA.*indexWeightsA + syncedCircSeqB.*syncedIndexWeightsB;
    indexWeightsAB = indexWeightsA + syncedIndexWeightsB;
    bitmaskAB = logical(indexWeightsAB);
    circSeqAB(bitmaskAB) = circSeqAB(bitmaskAB)./indexWeightsAB(bitmaskAB);
end

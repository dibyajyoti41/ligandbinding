function [barcode, bitmask, stdErrOfTheMean, indexWeights, clusterResultStruct, clusterBarcodeStruct] = extract_cluster_deliverables(consensusStruct, clusterKey)
    clusterKeys = consensusStruct.clusterKeys;
    clusterNumber = find(cellfun(@(s) strcmp(s, clusterKey), clusterKeys), 1, 'first');
    if isempty(clusterNumber)
        error(['Cluster key not found: ''', clusterKey, '''']);
    end
    clusterBarcodeStruct = consensusStruct.barcodeStructsMap(clusterKey);
    clusterResultStruct = consensusStruct.clusterResultStructs{clusterNumber};
    matAlignedBarcodes = cat(1, clusterResultStruct.alignedBarcodes{:});
    matAlignedBitmasks = cat(1, clusterResultStruct.alignedBarcodeBitmasks{:});
    matAlignedBarcodes(~matAlignedBitmasks) = NaN;
    barcode = clusterBarcodeStruct.barcode;
    indexWeights = clusterBarcodeStruct.indexWeights;
    bitmask = logical(indexWeights);
    barcodeLen = length(barcode);
    stdErrOfTheMean = NaN(1, barcodeLen);
    for idx=1:barcodeLen
        if bitmask(idx)
            vals = matAlignedBarcodes(:, idx);
            vals = vals(matAlignedBitmasks(:, idx));
            stdErrOfTheMean(idx) = std(vals)/sqrt(length(vals));
        end
    end
end
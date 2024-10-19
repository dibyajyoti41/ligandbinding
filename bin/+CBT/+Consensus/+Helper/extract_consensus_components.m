function [barcodeKeys, barcodes, barcodeBitmasks, flipTFs, circShifts] = extract_consensus_components(consensusKey, barcodeStructsMap)
    import CBT.Consensus.Helper.unfold_consensus;
    unfoldedComponents = unfold_consensus(consensusKey, barcodeStructsMap);
    numComponents = size(unfoldedComponents, 1);
    barcodeKeys = cell(numComponents,1);
    flipTFs = false(numComponents, 1);
    circShifts = zeros(numComponents, 1);
    for componentNum=1:numComponents
        [barcodeKeys{componentNum}, flipTFs(componentNum), circShifts(componentNum)] = unfoldedComponents{componentNum,:};
    end
    barcodeStructs = cellfun(@(x) barcodeStructsMap(x), barcodeKeys, 'UniformOutput', false);
    barcodes = cellfun(@(x) x.barcode, barcodeStructs, 'UniformOutput', false);
    barcodeBitmasks = cellfun(@(x) logical(x.indexWeights), barcodeStructs, 'UniformOutput', false);
end
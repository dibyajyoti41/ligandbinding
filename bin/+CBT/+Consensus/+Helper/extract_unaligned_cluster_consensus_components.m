function [unalignedBarcodes, unalignedBitmasks, clusterBarcodeKeys, clusterBarcodeAliases, otherBarcodeData, clusterBarcodeNums] = extract_unaligned_cluster_consensus_components(clusterConsensusData)
    consensusStruct = clusterConsensusData.details.consensusStruct;
    clusterBarcodeKeys = clusterConsensusData.clusterResultStruct.barcodeKeys;
    clusterBarcodeNums = str2double(clusterBarcodeKeys);

    inputs = consensusStruct.inputs;
    clusterBarcodeAliases = inputs.barcodeAliases(clusterBarcodeNums);
    unalignedBarcodes = inputs.barcodes(clusterBarcodeNums);
    unalignedBitmasks = inputs.barcodeBitmasks(clusterBarcodeNums);
    otherBarcodeData = [];
    if isfield(inputs, 'otherBarcodeData')
        inputs.otherBarcodeData(clusterBarcodeNums);
    end
end
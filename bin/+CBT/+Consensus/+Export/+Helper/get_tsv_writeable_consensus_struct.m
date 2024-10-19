function [dataStruct, columnFields, columnNames] = get_tsv_writeable_consensus_struct(clusterConsensusData)
    import CBT.Consensus.Helper.extract_aligned_cluster_consensus_components;
    [alignedBarcodes, alignedBitmasks, clusterBarcodeKeys, clusterBarcodeAliases] = extract_aligned_cluster_consensus_components(clusterConsensusData);
    barcodeColumnNames = strcat({'Aligned Barcode '}, clusterBarcodeKeys, ' (', clusterBarcodeAliases, ')');
    bitmaskColumnNames = strcat({'Aligned Bitmask '}, clusterBarcodeKeys, ' (', clusterBarcodeAliases, ')');
    numBarcodeColumns = size(barcodeColumnNames, 1);
    columnFields = cell(3 + numBarcodeColumns*2, 1);
    columnNames = cell(3 + numBarcodeColumns*2, 1);
    columnFields(1:3) = {'consensusBarcode'; 'consensusBitmask'; 'consensusIndexWeights'};
    columnNames(1:3) = strcat({'Cluster Consensus Barcode '; 'Cluster Consensus Bitmask '; 'Cluster Consensus Index Weights '}, clusterConsensusData.clusterKey);
    dataStruct = struct;
    dataStruct.consensusBarcode = clusterConsensusData.barcode;
    dataStruct.consensusBitmask = double(clusterConsensusData.bitmask);
    dataStruct.consensusIndexWeights = clusterConsensusData.indexWeights;
    for barcodeColumnNum=1:numBarcodeColumns
        barcodeFieldNum = (3 + (barcodeColumnNum - 1)*2) + 1;
        bitmaskFieldNum = barcodeFieldNum + 1;
        columnFields([barcodeFieldNum, bitmaskFieldNum]) = {...
            sprintf('Barcode_%s', clusterBarcodeKeys{barcodeColumnNum});...
            sprintf('Bitmask_%s', clusterBarcodeKeys{barcodeColumnNum})...
        };
        columnNames([barcodeFieldNum, bitmaskFieldNum]) = {...
            barcodeColumnNames{barcodeColumnNum};...
            bitmaskColumnNames{barcodeColumnNum}...
        };
        dataStruct.(columnFields{barcodeFieldNum}) = alignedBarcodes{barcodeColumnNum};
        dataStruct.(columnFields{bitmaskFieldNum}) = double(alignedBitmasks{barcodeColumnNum});
    end
end
function [tiffFilepath] = export_cluster_tiff(clusterConsensusData, clusterKey)
    import CBT.Consensus.Export.Helper.get_default_consensus_dirpath;
    
    timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
    defaultConsensusDirpath = get_default_consensus_dirpath();
    defaultTiffFilename = sprintf('AC_%s_%s.tif', strrep(strrep(clusterKey, '[', '('), ']', ')'), timestamp);
    defaultTiffFilepath = fullfile(defaultConsensusDirpath, defaultTiffFilename);
    [tiffFilename, tiffDirpath, ~] = uiputfile('*.tif', 'Save Automatic Consensus Cluster Image As', defaultTiffFilepath);
    if isequal(tiffDirpath, 0)
        return;
    end
    tiffFilepath = fullfile(tiffDirpath, tiffFilename);
    
    import CBT.Consensus.Helper.extract_aligned_cluster_consensus_components;
    [alignedBarcodes, alignedBitmasks] = extract_aligned_cluster_consensus_components(clusterConsensusData);
    consensusBarcode = clusterConsensusData.barcode;
    consensusBitmask = clusterConsensusData.bitmask;
    
    import CBT.Consensus.Core.gen_agg_consensus_barcodes_image;
    clusterConsensusImage = gen_agg_consensus_barcodes_image(alignedBarcodes, alignedBitmasks, consensusBarcode, consensusBitmask);
    imwrite(clusterConsensusImage, tiffFilepath);
end
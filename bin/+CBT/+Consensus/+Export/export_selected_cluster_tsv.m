function [] = export_selected_cluster_tsv(clusterConsensusData, clusterKey)
    import CBT.Consensus.Export.Helper.get_default_consensus_dirpath;
    import CBT.Consensus.Export.Helper.get_tsv_writeable_consensus_struct;
    import Fancy.IO.TSV.write_tsv;

    timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
    defaultConsensusDirpath = get_default_consensus_dirpath();
    defaultTsvFilename = sprintf('AC_%s_%s.tsv', strrep(strrep(clusterKey, '[', '('), ']', ')'), timestamp);
    defaultTsvFilepath = fullfile(defaultConsensusDirpath, defaultTsvFilename);
    [tsvFilename, tsvDirpath] = uiputfile('*.tsv', 'Save Automatic Consensus Cluster TSV Spreadsheet As', defaultTsvFilepath);
    if isequal(tsvDirpath, 0)
        return;
    end
    tsvFilepath = fullfile(tsvDirpath, tsvFilename);
    [dataStruct, columnFields, columnNames] = get_tsv_writeable_consensus_struct(clusterConsensusData);
    write_tsv(tsvFilepath, dataStruct, columnFields, columnNames);
end
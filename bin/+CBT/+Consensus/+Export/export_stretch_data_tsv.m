function [] = export_stretch_data_tsv(clusterConsensusData, clusterKey)
    import CBT.Consensus.Export.Helper.get_default_consensus_dirpath;
    import CBT.Consensus.Export.Helper.get_tsv_writeable_consensus_stretching_struct;
    import Fancy.IO.TSV.write_tsv;
    
    timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
    defaultConsensusDirpath = get_default_consensus_dirpath();
    defaultConsensusTsvFilename = sprintf('AC_stretch_%s_%s.tsv', strrep(strrep(clusterKey, '[', '('), ']', ')'), timestamp);
    defaultConsensusTsvFilepath = fullfile(defaultConsensusDirpath, defaultConsensusTsvFilename);
    [consensusTsvFilename, consensusTsvDirpath] = uiputfile('*.tsv', 'Save Cluster Stretch Data TSV Spreadsheet As', defaultConsensusTsvFilepath);
    if isequal(consensusTsvDirpath, 0)
        return;
    end
    consensusTsvFilepath = fullfile(consensusTsvDirpath, consensusTsvFilename);
    [dataStruct, columnFields, columnNames] = get_tsv_writeable_consensus_stretching_struct(clusterConsensusData);
    write_tsv(consensusTsvFilepath, dataStruct, columnFields, columnNames);
    
    fprintf('Saved cluster stretch data spreadsheet for cluster ''%s'' to ''%s''\n', clusterKey, consensusTsvFilepath);

end
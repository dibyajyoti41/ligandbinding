function [matFilepath] = export_cluster_mat(clusterConsensusData, clusterKey) %#ok<INUSL>
    timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
    import CBT.Consensus.Export.Helper.get_default_consensus_dirpath;
    defaultConsensusDirpath = get_default_consensus_dirpath();
    defaultMatFilename = sprintf('AC_%s_%s.mat', strrep(strrep(clusterKey, '[', '('), ']', ')'), timestamp);
    defaultMatFilepath = fullfile(defaultConsensusDirpath, defaultMatFilename);
    [matFilename, matDirpath] = uiputfile('*.mat', 'Save Automatic Consensus Cluster Data As', defaultMatFilepath);
    if isequal(matDirpath, 0)
        return;
    end
    matFilepath = fullfile(matDirpath, matFilename);
    save(matFilepath, 'clusterConsensusData');
end

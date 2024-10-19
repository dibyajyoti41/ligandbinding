function [clusterConsensusData, consensusFilepath] = load_cluster_consensus(consensusFilepath)
    if nargin < 1
        consensusFilepath = [];
    end
    clusterConsensusData = [];
    import OptMap.DataImport.try_prompt_consensus_filepaths;
    [~, consensusFilepaths] = try_prompt_consensus_filepaths([], false);
    if isempty(consensusFilepaths)
        return;
    end
    consensusFilepath = consensusFilepaths{1};
    if isempty(consensusFilepath)
        return;
    end
    consensusStruct = load(consensusFilepath, 'clusterConsensusData');
    if isfield(consensusStruct, 'clusterConsensusData')
        clusterConsensusData = consensusStruct.clusterConsensusData;
    end
end
function [] = load_consensus_results(ts, consensusStruct)
    if nargin < 2
        import CBT.Consensus.Import.load_cluster_consensus;
        clusterConsensusData = load_cluster_consensus();
        if isempty(clusterConsensusData)
            return;
        end
        consensusStruct = clusterConsensusData.details.consensusStruct;
    end
    if isempty(consensusStruct)
        return;
    end

    assignin('base', 'consensusStruct', consensusStruct);
    import CBT.Consensus.UI.generate_consensus_ui;
    generate_consensus_ui(ts, consensusStruct);
end
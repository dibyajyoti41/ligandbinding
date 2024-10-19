function [] = generate_consensus_ui(ts, consensusStruct)
    import CBT.Consensus.Helper.make_cluster_structs;
    [clusterKeys, clusterConsensusDataStructs] = make_cluster_structs(consensusStruct);

    import CBT.Consensus.UI.launch_export_ui;
    launch_export_ui(ts, clusterKeys, clusterConsensusDataStructs)

    import CBT.Consensus.UI.plot_clusters_concentrically;
    hTabConcentric = plot_clusters_concentrically(ts, clusterKeys, clusterConsensusDataStructs);

    import CBT.Consensus.UI.plot_clusters_linearly;
    plot_clusters_linearly(ts, consensusStruct);

    hTabPairwiseConsensusHistory = ts.create_tab('Pairwise Consensus History');
    ts.select_tab(hTabPairwiseConsensusHistory);
    hPanelPairwiseConsensusHistory = uipanel(hTabPairwiseConsensusHistory);
    import CBT.Consensus.UI.plot_pairwise_consensus_history;
    plot_pairwise_consensus_history(consensusStruct, hPanelPairwiseConsensusHistory);

    hTabConsensusDendros = ts.create_tab('Consensus Dendrograms');
    ts.select_tab(hTabConsensusDendros);
    hPanelConsensusDendros = uipanel(hTabConsensusDendros);
    import CBT.Consensus.UI.plot_consensus_dendrograms;
    plot_consensus_dendrograms(consensusStruct, hPanelConsensusDendros);

    ts.select_tab(hTabConcentric);
end
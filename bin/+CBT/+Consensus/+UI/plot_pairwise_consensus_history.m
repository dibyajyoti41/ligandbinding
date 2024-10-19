function [] = plot_pairwise_consensus_history(consensusStruct, hPanel)
    hPanel = uipanel('Parent', hPanel);

    import Fancy.UI.FancyTabs.TabbedScreen;
    ts = TabbedScreen(hPanel);
    
    import CBT.Consensus.UI.plot_pairwise_consensus_history_helper;
    plot_pairwise_consensus_history_helper(ts, consensusStruct.finalConsensusKey, consensusStruct.barcodeStructsMap);
end
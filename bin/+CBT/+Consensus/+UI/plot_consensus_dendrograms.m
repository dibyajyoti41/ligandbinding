function [] = plot_consensus_dendrograms(consensusStruct, hPanel)
    import Fancy.UI.FancyTabs.TabbedScreen;
    ts = TabbedScreen(hPanel);

    barcodeAliases = consensusStruct.inputs.barcodeAliases;
    numBarcodes = length(consensusStruct.inputs.barcodes);
    commonLength = length(consensusStruct.inputs.barcodes{1});
    if not(isfield(consensusStruct.inputs, 'clusterThresholdScore'))
        warning('Old consensus information -- needs to be fixed with OptMap.Scripts.Compatibility.fix_consensus_files');
        return;
    end
    clusterThresholdScore = consensusStruct.inputs.clusterThresholdScore;
    bestPossibleScore = sqrt(commonLength);
    clusterThresholdScoreNormalized = clusterThresholdScore/bestPossibleScore;
    consensusMergingTree = consensusStruct.consensusMergingTree;
    clusterAssignmentsMatrix = consensusStruct.clusterAssignmentsMatrix;

    
    dendroDefaultOrderTitleStr = 'Consensus Dendrogram (Default Order)';
    hDendroDefaultOrderTab = ts.create_tab(dendroDefaultOrderTitleStr);
    ts.select_tab(hDendroDefaultOrderTab);
    
    import CBT.Consensus.UI.plot_dendrogram;
    plot_dendrogram(hDendroDefaultOrderTab, numBarcodes, consensusMergingTree, clusterThresholdScoreNormalized, [], barcodeAliases, commonLength);

    
    dendroReorderedTitleStr = 'Consensus Dendrogram (Reordered)';
    hDendroReorderedTab = ts.create_tab(dendroReorderedTitleStr);
    ts.select_tab(hDendroReorderedTab);
    [~, leafReordering] = sort(clusterAssignmentsMatrix(:));
    leafReordering = leafReordering';

    import CBT.Consensus.UI.plot_dendrogram;
    plot_dendrogram(hDendroReorderedTab, numBarcodes, consensusMergingTree, clusterThresholdScoreNormalized,  leafReordering, barcodeAliases, commonLength);
end
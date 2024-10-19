function [hTabConcentric, tsConcentric] = plot_clusters_concentrically(ts, clusterKeys, clusterConsensusDataStructs)

    import Fancy.UI.FancyTabs.TabbedScreen;

    hTabConcentric = ts.create_tab('Concentric Plots');
    hPanelConcentric = uipanel('Parent', hTabConcentric);
    tsConcentric = TabbedScreen(hPanelConcentric);

    
    numBarcodesInClusters = cellfun(@(clusterKey) sum(clusterKey == ','), clusterKeys) + 1;
    [~, sortOrder] = sort(numBarcodesInClusters, 'descend');
    numClusters = length(clusterKeys);
    
    import CBT.Consensus.UI.plot_cluster_concentrically;
    for clusterNum = 1:numClusters
        clusterIdx = sortOrder(clusterNum);
        clusterKey = clusterKeys{clusterIdx};
        clusterConsensusData = clusterConsensusDataStructs{clusterIdx};
        hTabCluster = tsConcentric.create_tab(sprintf('%s', clusterKey));
        hPanelCluster = uipanel('Parent', hTabCluster);
        tsCluster = TabbedScreen(hPanelCluster);
        plot_cluster_concentrically(tsCluster, clusterConsensusData);
    end
end
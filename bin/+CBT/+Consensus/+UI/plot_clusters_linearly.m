function [hTabLinear] = plot_clusters_linearly(ts, consensusStruct)
    hTabLinear = ts.create_tab('Linear Plots');
    hPanelLinear = uipanel(hTabLinear);

    import Fancy.UI.FancyTabs.TabbedScreen;
    tsLinear = TabbedScreen(hPanelLinear);

    clusterKeys = consensusStruct.clusterKeys;
    barcodeStructsMap = consensusStruct.barcodeStructsMap;
    clusterResultStructs = consensusStruct.clusterResultStructs;

    numClusters = size(clusterKeys, 1);
    if numClusters == 0
        return;
    end
    import CBT.Consensus.UI.plot_cluster_linearly;
    for clusterNum = 1:numClusters
        clusterKey = clusterKeys{clusterNum};
        clusterStruct = barcodeStructsMap(clusterKey);
        clusterResultStruct = clusterResultStructs{clusterNum};
        hTabLinearCluster = tsLinear.create_tab(clusterKey);
        hAxis = axes('Units', 'normal', 'Position', [0 0 1 1], 'Parent', hTabLinearCluster);
        plot_cluster_linearly(hAxis, clusterKey, clusterStruct, clusterResultStruct, barcodeStructsMap);
    end
end
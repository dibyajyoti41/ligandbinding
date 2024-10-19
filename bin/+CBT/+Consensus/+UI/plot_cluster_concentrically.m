function [] = plot_cluster_concentrically(ts, clusterConsensusData, plotUnalignedBarcodesTF)
    if nargin < 3
        plotUnalignedBarcodesTF = false;
    end
    clusterKey = clusterConsensusData.clusterKey;
    clusterName = clusterKey;

    
    if plotUnalignedBarcodesTF
        import CBT.Consensus.Helper.extract_unaligned_cluster_consensus_components;
        [unalignedBarcodes, unalignedBitmasks, clusterBarcodeKeys, clusterBarcodeAliases] = extract_unaligned_cluster_consensus_components(clusterConsensusData);
        clusterBarcodeNames1 = strcat({'Unaligned Barcode '}, clusterBarcodeKeys, {' ('}, clusterBarcodeAliases, {')'});

        import Barcoding.Helpers.sanitize_barcodes;
        sanitizedBarcodesUnaligned = sanitize_barcodes(unalignedBarcodes, unalignedBitmasks);
        tabTitleUnaligned = 'Unaligned';
        titleUnaligned = {'Unaligned Barcodes in Consensus Cluster (Concentric Plot)', clusterName};
        hTabUnaligned = ts.create_tab(tabTitleUnaligned);
        ts.select_tab(hTabUnaligned);
        hPanelUnaligned = uipanel('Parent', hTabUnaligned);
        hAxisUnalignes = axes(...
            'Units', 'normal', ...
            'Position', [0, 0.1, 0.6, 0.8], ...
            'Parent', hPanelUnaligned ...
            );

        uitable(hPanelUnaligned,...
            'Data', clusterBarcodeNames1,...
            'ColumnName', {'Barcode Names (Numbered from center outwards)'},...
            'ColumnWidth', {400},...
            'Units', 'normal', 'Position', [0.6, 0, 0.4, 1]);
        import Barcoding.Visualizing.plot_circular_barcodes_concentrically;
        plot_circular_barcodes_concentrically(hAxisUnalignes, [sanitizedBarcodesUnaligned; {nan(size(sanitizedBarcodesUnaligned{1}))}]);
        title(titleUnaligned);
    end



    import CBT.Consensus.Helper.extract_aligned_cluster_consensus_components;
    [alignedBarcodes, alignedBitmasks, clusterBarcodeKeys, clusterBarcodeAliases] = extract_aligned_cluster_consensus_components(clusterConsensusData);
    clusterBarcodeNamesAligned = strcat({'Aligned Barcode '}, clusterBarcodeKeys, {' ('}, clusterBarcodeAliases, {')'});

    import Barcoding.Helpers.sanitize_barcodes;
    sanitizedBarcodesAligned = sanitize_barcodes(alignedBarcodes, alignedBitmasks);
    tabTitleAligned = 'Aligned';
    titleAligned = {'Aligned Barcodes in Consensus Cluster (Concentric Plot)', clusterName};



    hTabAligned = ts.create_tab(tabTitleAligned);
    ts.select_tab(hTabAligned);
    hPanelAligned = uipanel('Parent', hTabAligned);
    hAxisAligned = axes(...
        'Units', 'normal', ...
        'Position', [0, 0.1, 0.6, 0.8], ...
        'Parent', hPanelAligned ...
        );

    uitable(hPanelAligned,...
        'Data', clusterBarcodeNamesAligned,...
        'ColumnName', {'Barcode Names (Numbered from center outwards)'},...
        'ColumnWidth', {400},...
        'Units', 'normal', 'Position', [0.6, 0, 0.4, 1]);
    sanitizedBarcodes = [sanitizedBarcodesAligned; {sanitize_barcodes(clusterConsensusData.barcode, clusterConsensusData.bitmask)}];
    sanitizedBarcodes = cellfun(@(v) (v - nanmean(v(:)))./nanstd(v(:)), sanitizedBarcodes, 'UniformOutput', false);
    import Barcoding.Visualizing.plot_circular_barcodes_concentrically;
    plot_circular_barcodes_concentrically(hAxisAligned, sanitizedBarcodes);
    title(hAxisAligned, titleAligned);
end
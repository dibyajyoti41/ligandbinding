function [] = plot_pairwise_consensus_history_helper(ts, consensusStructKey, barcodeStructsMap)
    import Barcoding.Reorienting.reorient_barcode_and_bitmask;

    consensusStruct = barcodeStructsMap(consensusStructKey);
    if isempty(consensusStruct.alias)
        consensusLabel = consensusStructKey;
    else
        consensusLabel = [consensusStructKey, ' (',  strrep(consensusStruct.alias, '_', '\_'), ')'];
    end
    keys = cell(2,1);
    flipTFs = zeros(2,1);
    circShifts = zeros(2,1);
    barcodes = cell(2,1);
    barcodeBitmasks = cell(2,1);
    reorientedBarcodes = cell(2,1);
    reorientedBarcodeBitmasks = cell(2,1);
    parents = consensusStruct.parents;
    if isempty(parents)
        return;
    end
    hTab = ts.create_tab(consensusStructKey);
    ts.select_tab(hTab);
    hAxis = axes('Units','normal', 'Position', [0 0 1 1], 'Parent', hTab);
    selfBarcode = consensusStruct.barcode;
    selfBitmask = logical(consensusStruct.indexWeights);
    values = selfBarcode;
    values(~selfBitmask) = NaN;
    consensusPlot = plot(hAxis, values, '-o');
    hold(hAxis, 'on');
    plots = gobjects(2,1);
    labels = cell(2,1);
    for parentComponentNum=1:2
        parent = parents{parentComponentNum};
        [keys{parentComponentNum}, flipTFs(parentComponentNum), circShifts(parentComponentNum)] = parent{:};
        barcodeStruct = barcodeStructsMap(keys{parentComponentNum});
        if isempty(barcodeStruct.alias)
            labels{parentComponentNum} = keys{parentComponentNum};
        else
            labels{parentComponentNum} = [keys{parentComponentNum}, ' (', strrep(barcodeStruct.alias, '_', '\_'), ')'];
        end
        barcodes{parentComponentNum} = barcodeStruct.barcode;
        barcodeBitmasks{parentComponentNum} = logical(barcodeStruct.indexWeights);
        [reorientedBarcodes{parentComponentNum}, reorientedBarcodeBitmasks{parentComponentNum}] = reorient_barcode_and_bitmask(...
            barcodes{parentComponentNum}, ...
            barcodeBitmasks{parentComponentNum}, ...
            flipTFs(parentComponentNum), ...
            circShifts(parentComponentNum));
        values = reorientedBarcodes{parentComponentNum};
        values(~reorientedBarcodeBitmasks{parentComponentNum}) = NaN;
        plots(parentComponentNum) = plot(hAxis, values);
    end
    legend(hAxis, [consensusPlot; plots], [consensusLabel; labels]);
    
    import CBT.Consensus.UI.plot_pairwise_consensus_history_helper;
    for parentComponentNum=1:2
        plot_pairwise_consensus_history_helper(ts, keys{parentComponentNum}, barcodeStructsMap);
    end
    ts.select_tab(hTab);
end
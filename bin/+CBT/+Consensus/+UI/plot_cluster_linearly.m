function [] = plot_cluster_linearly(hAxis, clusterKey, clusterStruct, clusterResultStruct, barcodeStructsMap)
    numComponents = size(clusterResultStruct.barcodeKeys, 1);
    hPlots = gobjects(numComponents, 1);

    selfBarcode = clusterStruct.barcode;
    selfBitmask = logical(clusterStruct.indexWeights);
    if isempty(clusterStruct.alias)
        clusterLabel = clusterKey;
    else
        clusterLabel = [clusterKey, ' (', strrep(clusterStruct.alias, '_', '\_'), ')'];
    end
    values = selfBarcode;
    values(~selfBitmask) = NaN;
    values = (values - nanmean(values(:)))./nanstd(values(:));
    hPlotSelf = plot(hAxis, values, '-o');
    hold(hAxis, 'on');
    labels = cell(numComponents, 1);
    for componentNum = 1:numComponents
        key = clusterResultStruct.barcodeKeys{componentNum};
        componentStruct = barcodeStructsMap(key);
        if isempty(componentStruct.alias)
            labels{componentNum} = key;
        else
            labels{componentNum} = [key, ' (', strrep(componentStruct.alias, '_', '\_'), ')'];
        end
        values = clusterResultStruct.alignedBarcodes{componentNum};
        values(~clusterResultStruct.alignedBarcodeBitmasks{componentNum}) = NaN;
        values = (values - nanmean(values(:)))./nanstd(values(:));
        hPlots(componentNum) = plot(hAxis, values);
    end
    legend(hAxis, [hPlotSelf; hPlots], [clusterLabel; labels]);
end
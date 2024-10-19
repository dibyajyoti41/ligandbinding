function [] = plot_dendrogram(hParent, numBarcodes, consensusMergingTreeMatrix, clusterThresholdScoreNormalized, leafReordering, barcodeAliases, commonLength)

    if nargin < 5
        leafReordering = [];
    end
    if nargin < 6
        barcodeAliases = cell(numBarcodes, 1);
    end
    
    hAxis = axes(...
        'Parent', hParent, ...
        'Units', 'normal', ...
        'Position', [.1 .1 .6 .8], ...
        'FontSize', 15);
    if not(isempty(leafReordering))
        hDendroLines = dendrogram(consensusMergingTreeMatrix, 0, 'ColorThreshold', (1 - clusterThresholdScoreNormalized), 'Reorder', leafReordering, 'CheckCrossing', false);
    else
        hDendroLines = dendrogram(consensusMergingTreeMatrix, 0, 'ColorThreshold', (1 - clusterThresholdScoreNormalized));
    end
    set(hDendroLines, 'LineWidth', 2);
    y = (1 - clusterThresholdScoreNormalized)*[1 1];
    hold(hAxis, 'on');
    plot(hAxis, xlim(hAxis), y, 'b:');
    xlabel(hAxis, 'Barcode Number', 'FontSize', 14);

    xTickLabels = get(hAxis, 'XTickLabel');
    if ischar(xTickLabels)
        xTickLabels = mat2cell(xTickLabels, ones(size(xTickLabels, 1), 1), size(xTickLabels, 2));
    end
    
    yTickLabels = get(hAxis, 'YTickLabel');
    if ischar(yTickLabels)
        yTickLabels = mat2cell(yTickLabels, ones(size(yTickLabels, 1), 1), size(yTickLabels, 2));
    end
    
    yTickLabels = cellfun(@(y) num2str(1 - str2double(y)), yTickLabels, 'UniformOutput', false);
    xTickLabelsNumeric = cellfun(@str2double, xTickLabels);
    xTickLabels = arrayfun(@(x) num2str(x), xTickLabelsNumeric, 'UniformOutput', false);
    set(hAxis, 'FontSize', 12);
    set(hAxis, 'XTickLabel', xTickLabels);
    set(hAxis, 'YTickLabel', yTickLabels);
    set(hAxis, 'XTickLabelRotation', -90);
    ylabel(hAxis, ['$\frac{1}{\sqrt{N}} S_{ab}(\hat{f}, \hat{d})$ with $N = ', num2str(commonLength), '$'],'Interpreter','latex', 'FontSize', 18);

    hold(hAxis, 'off');

    tableColumnNames = {'Barcode Aliases'};
    tableRowNames = xTickLabelsNumeric;
    tableColumnWidths = {400};
    tableData = barcodeAliases(xTickLabelsNumeric);
    uitable(hParent,...
        'Data', tableData,...
        'ColumnName', tableColumnNames,...
        'RowName', tableRowNames,...
        'ColumnWidth', tableColumnWidths,...
        'Units', 'normal', 'Position', [.75 .1 .2 .8]);
end
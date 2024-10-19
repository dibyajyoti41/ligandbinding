function [] = position_ui_elems_in_grid(hControls, maxNumControlItemsPerRow, paddingNrm)
    import Fancy.Utils.merge_structs;
    import Fancy.UI.FancyPositioning.FancyGrid.generate_grid;
    import Fancy.UI.FancyPositioning.FancyGrid.get_grid_positions_nrm;
    import Fancy.UI.FancyPositioning.reposition_nrm;
    
    if nargin < 3
        paddingNrm = struct();
    end
    defaultPaddingNrm = struct('top', 0.0, 'bottom', 0.0, 'left', 0.0, 'right', 0.0);
    paddingNrm = merge_structs(defaultPaddingNrm, paddingNrm);
    
    % create grid labels and coordinates
    numControlItems = numel(hControls);
    numPadlessGridCols = min(numControlItems, maxNumControlItemsPerRow);
    numPadlessGridRows = ceil(double(numControlItems)/double(numPadlessGridCols));
    
    [padlessGridLabels, rowBoundaryPosFromTopNrm, colBoundaryPosFromLeftNrm] = generate_grid(numPadlessGridRows, numPadlessGridCols);
    padlessGridLabels(padlessGridLabels > numControlItems) = 0;

    % add padding to grid labels and coordinates
    
    paddedGridLabels = round(imresize(padlessGridLabels, 3*size(padlessGridLabels), 'nearest'));
    paddingMat = ordfilt2(paddedGridLabels, 1, true(3)) ~= ordfilt2(paddedGridLabels, 9, true(3));
    paddedGridLabels(paddingMat) = 0;
    rowBoundaryPosFromTopNrm = cumsum([rowBoundaryPosFromTopNrm(1), feval(@(tmp) [tmp{:}], arrayfun(@(x) [paddingNrm.top, x - paddingNrm.top - paddingNrm.bottom, paddingNrm.bottom], diff(rowBoundaryPosFromTopNrm), 'UniformOutput', false))]);
    colBoundaryPosFromLeftNrm = cumsum([colBoundaryPosFromLeftNrm(1), feval(@(tmp) [tmp{:}], arrayfun(@(x) [paddingNrm.left, x - paddingNrm.left - paddingNrm.right, paddingNrm.right], diff(colBoundaryPosFromLeftNrm), 'UniformOutput', false))]);
    numTmpRowsNew = size(paddedGridLabels, 1);
    for newTmpRowNum = numTmpRowsNew:-1:2
        if isequal(paddedGridLabels(newTmpRowNum, :), paddedGridLabels(newTmpRowNum - 1, :))
            paddedGridLabels = paddedGridLabels([1:(newTmpRowNum - 1), (newTmpRowNum + 1):end], :);
            rowBoundaryPosFromTopNrm(newTmpRowNum) = [];
        end
    end

    numTmpColsNew = size(paddedGridLabels, 2);
    for newTmpColNum = numTmpColsNew:-1:2
        if isequal(paddedGridLabels(:, newTmpColNum), paddedGridLabels(:, newTmpColNum - 1))
            paddedGridLabels = paddedGridLabels(:, [1:(newTmpColNum - 1), (newTmpColNum  + 1):end]);
            colBoundaryPosFromLeftNrm(newTmpColNum) = [];
        end
    end
    
    % reposition controls in accordance with labels and coordinates
    [cellVectChildPosArrNrm] = get_grid_positions_nrm(paddedGridLabels, rowBoundaryPosFromTopNrm, colBoundaryPosFromLeftNrm);
    reposition_nrm(hControls, cellVectChildPosArrNrm);
end
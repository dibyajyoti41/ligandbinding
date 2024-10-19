function [cellVectChildPosArrNrm] = get_grid_positions_nrm(gridLabels, rowBoundaryPosFromTopNrm, colBoundaryPosFromLeftNrm)
    validateattributes(gridLabels, {'numeric'}, {'nonnegative', 'integer', '2d'}, 1);
    gridSz = double(size(gridLabels));
    numGridRows = gridSz(1);
    numGridCols = gridSz(2);
    if nargin < 2
        rowBoundaryPosFromTopNrm = linspace(0, 1, numGridRows + 1);
    else
        validateattributes(rowBoundaryPosFromTopNrm, {'numeric'}, {'vector', 'numel', numGridRows + 1}, 2);
    end
    if nargin < 3
        colBoundaryPosFromLeftNrm = linspace(0, 1, numGridCols + 1);
    else
        validateattributes(colBoundaryPosFromLeftNrm, {'numeric'}, {'vector', 'numel', numGridCols + 1}, 3);
    end

    % adjust for coordinates increasing from bottom (0) to top (1)
    %  so they align with rows in a mtrix that go from top  to
    %    bottom
    numLabels = max(gridLabels(:));
    gridLabelsUD = flipud(gridLabels);
    rowBoundaryPosFromBottomNrm = flipud(1 - rowBoundaryPosFromTopNrm(:));
    falseGrid = false(gridSz);
    panelLeftPosNrm = NaN(numLabels, 1);
    panelBottomPosNrm = NaN(numLabels, 1);
    panelWidthsNrm = NaN(numLabels, 1);
    panelHeightsNrm = NaN(numLabels, 1);
    labelsFound = false(numLabels, 1);
    for labelNum = 1:numLabels
        currLabel = (gridLabelsUD == labelNum);
        if any(currLabel(:))
            labelsFound(labelNum) = true;
            [rowStart, colStart] = ind2sub(gridSz, find(currLabel, 1, 'first'));
            [rowEnd, colEnd] = ind2sub(gridSz, find(currLabel, 1, 'last'));
            currLabelValid = falseGrid;
            currLabelValid(rowStart:rowEnd, colStart:colEnd) = true;
            if not(isequal(currLabel, currLabelValid))
                error('Invalid gridLabel matrix: same-labeled areas must be rectangular and contiguous');
            end
            panelBottomPosNrm(labelNum) = rowBoundaryPosFromBottomNrm(rowStart);
            panelLeftPosNrm(labelNum) = colBoundaryPosFromLeftNrm(colStart);
            panelHeightsNrm(labelNum) = rowBoundaryPosFromBottomNrm(rowEnd + 1) - panelBottomPosNrm(labelNum);
            panelWidthsNrm(labelNum) = colBoundaryPosFromLeftNrm(colEnd + 1) - panelLeftPosNrm(labelNum);
        end
    end
    cellVectChildPosArrNrm = [panelLeftPosNrm, panelBottomPosNrm, panelWidthsNrm, panelHeightsNrm];
    cellVectChildPosArrNrm = mat2cell(cellVectChildPosArrNrm, ones([1, size(cellVectChildPosArrNrm, 1)]), size(cellVectChildPosArrNrm, 2));
    cellVectChildPosArrNrm(~labelsFound) = cell(sum(~labelsFound), 1);
end
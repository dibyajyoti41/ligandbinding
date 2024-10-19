function [gridLabels, rowBoundaryPosFromTopNrm, colBoundaryPosFromLeftNrm] = generate_grid(numGridRows, numGridCols)
    gridLabels = permute(reshape(1:(numGridRows * numGridCols), [numGridCols, numGridRows]), [2, 1]);
    rowBoundaryPosFromTopNrm = linspace(0, 1, numGridRows + 1);
    colBoundaryPosFromLeftNrm = linspace(0, 1, numGridCols + 1);
end
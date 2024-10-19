function pathPlotMat = gen_labeled_paths_mat(pathsArr, matSize)  
    % GEN_LABELED_PATHS_MAT - Maps out the paths found in a kymograph.
    %   Creates a sparse array of the same size as the kymograph in which 
    %   the features were found, containing non-zero elements where pixels 
    %   belong to a path. The value of these non-zero elements are integer 
    %   values corresponing to the idex of the path.
    %   The results of this function can be drawn directly using 
    %   imshow(pathPlotMat)
    %
    % Inputs:
    %  paths
    %   An array, where each column represents one path, and contains all 
    %   pixels belonging to the path.
    %  matSize
    %   the size of the kymograph in which the features were found
    %
    % Outputs:
    %  pathArr
    %   Sparse array of the same size as the kymograph in which the 
    %   features were found, where non-zero elements correspond to pixels
    %   being crossed by path.
    %
    % 
    % Authors:
    %  Henrik Nordanger

    numRows = size(pathsArr, 1);
    numPaths = size(pathsArr, 2);

    pathPlotMat = zeros(matSize);

    for pathNum = 1:numPaths
        for rowIdx = 1:numRows
            colIdx = pathsArr(rowIdx, pathNum);
            pathPlotMat(rowIdx, colIdx) = pathNum;
        end
    end
end
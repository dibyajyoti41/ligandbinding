function [stretchFactorsMat] = compute_stretch_factors(pathsColIdxs, imgSz)
    % COMPUTE_STRETCH_FACTORS - computes the stretch factors for the
    %    image
    %
    % Inputs:
    %   pathsColIdxs
    %     the output column indices of the paths for each feature
    %      to be straightened in the stretched output
    %   imgSz
    %     the size of the image
    % 
    % Outputs:
    %   stretchFactorsMat
    %     matrix containing the horizontal stretch factor of each pixel
    % 
    % Authors:
    %   Saair Quaderi (2017-02): Separated from apply_horizontal_stretching
    %      code
    %   Henrik Nordanger (2016-10): Generalized to stretch more than
    %      one path
    %   Charleston Noble

    numRows = imgSz(1);
    numCols = imgSz(2);
    validateattributes(pathsColIdxs, {'numeric'}, {'positive', 'integer', 'nrows', numRows}, 1);

    %The paths are sorted.
    pathsColIdxs = sortrows(pathsColIdxs');
    pathsColIdxs = pathsColIdxs';

    %The number of features is obtained.
    numFeatures = size(pathsColIdxs, 2);

    %The mean position of each feature is determined.
    xMean = zeros(1, numFeatures);
    for featureNum = 1:numFeatures
        xMean(featureNum) = mean(pathsColIdxs(:, featureNum));
    end

    %Array for containing the stretch factor for each "gap" between
    %features is pre-allocated.
    regionalStretchFactors = zeros(numRows, numFeatures + 1);
    stretchFactorsMat = ones(numRows, numCols);

    for rowNum = 1:numRows
        stretchRow = pathsColIdxs(rowNum, :); 

        % The first stretch factor for the current row is determined
        regionalStretchFactors(rowNum, 1) = xMean(1) / stretchRow(1);
        stretchFactorsMat(rowNum, (1:stretchRow(1))) = regionalStretchFactors(rowNum, 1);

        % The stretch factors between each of the features are determined
        for featureNum = 2:numFeatures
            regionalStretchFactors(rowNum, featureNum) = (xMean(featureNum) - xMean(featureNum - 1)) / (stretchRow(featureNum) - stretchRow(featureNum - 1));
            stretchFactorsMat(rowNum, (stretchRow(featureNum-1) + 1):stretchRow(featureNum)) = regionalStretchFactors(rowNum, featureNum);    
        end

        % The last stretch factor for the current row is determined
        regionalStretchFactors(rowNum, end) = (numCols - xMean(end)) / (numCols - stretchRow(end));
        stretchFactorsMat(rowNum, (stretchRow(end) + 1):end) = regionalStretchFactors(rowNum, end);
    end
end
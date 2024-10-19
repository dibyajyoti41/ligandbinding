function foundFeatures = find_k_features(prefilteredImg, maxMovementPx, typicalFeatureWidth, maxNumFeaturesSoughtK)
    % FIND_K_FEATURES - detects a single "feature" in the WPAlign
    %	algorithm (see the WPAlign paper for details)
    %
    % Inputs:
    %  prefilteredImg
    %    the image to detect the feature in (pre-filtered to express
    %    features and to have values between -1 and 1)
    %  maxMovementPx
    %    the number of pixels to the left/right that the
    %      feature can jump in each step
    %  typicalFeatureHalfwidth
    %    the number of pixels to the left/right of a found feature to
    %    ignore
    %  maxNumFeaturesSoughtK
    %    the maximum number of features to find
    %
    %
    % Outputs:
    %   foundFeatures
    %     cells where
    %       entries in the first column contains the path of the feature
    %         (cell of the positions for a feature)
    %       entries in the second column are the path distances (costs)
    %         of the featur
    % 
    % Authors:
    %   Henrik Nordanger (2016-10): eliminated crossing paths problem
    %   Saair Quaderi (2015-11): refactoring
    %   Charleston Noble

    prefilteredImgNonnan = prefilteredImg;
    prefilteredImgNonnan(isnan(prefilteredImgNonnan)) = 0;
    validateattributes(prefilteredImgNonnan, {'numeric'}, {'<=', 1, '>=' -1}, 1);
    validateattributes(maxMovementPx, {'numeric'}, {'scalar', 'nonnan', 'real', 'nonnegative'}, 2);
    maxMovementPx = floor(maxMovementPx);

    imgSize = size(prefilteredImg);

    imgUp = inf(imgSize);
    imgDown = inf(imgSize);

    upMask = prefilteredImg > 0;
    downMask = prefilteredImg < 0;

    imgUp(upMask) = 1 - prefilteredImg(upMask);
    imgDown(downMask) = 1 + prefilteredImg(downMask);

    weightImgs = {imgUp; imgDown};

    numRows = imgSize(1);
    numCols = imgSize(2);

    idxs = zeros(1,maxMovementPx);
    idxs(1) = 1;

    P = maxMovementPx - 1;
    halfP = P/2;
    M = numRows * (numCols + P) + 2;    
    
    trimPx = (typicalFeatureWidth - 1) / 2;

    % Creating I with padding connections
    iPart1 = ones(1, numCols);    
    iPart2 = 1 + (1:((numRows - 1)*(numCols + P)));

    paddedImgColIdxsNegOne = 1:P;
    paddedImgColIdxsNegOne(paddedImgColIdxsNegOne > halfP) = paddedImgColIdxsNegOne(paddedImgColIdxsNegOne > halfP) + numCols;

    removeIdxsA = (numCols + P).*(repmat((1:(numRows - 1))', 1, P) - 1);
    removeIdxsB = repmat(paddedImgColIdxsNegOne, numRows - 1, 1);
    removeIdxs = removeIdxsA + removeIdxsB;
    iPart2(removeIdxs(:)) = [];    
    iPart2 = iPart2(cumsum(repmat(idxs, 1, length(iPart2))));
    iPart3 = (1:numCols) + (numRows - 1)*(numCols + P) + halfP + 1;
    I = [iPart1, iPart2, iPart3];

    % Creating J with padding connections
    jPart1 = (halfP + 1) + (1:numCols);
    x = 1:(numRows - 1);
    x = (numCols + P).*x;
    x = x + 2;
    xRep = repmat(x, numCols, 1);
    jumpMat = cumsum(ones(numCols, length(x))) - 1;
    idx = xRep + jumpMat;
    idx = idx(:)';
    idxRep = repmat(idx, maxMovementPx, 1);
    jumpMat = cumsum(ones(maxMovementPx, length(idxRep))) - 1;
    jPart2 = idxRep + jumpMat;
    jPart2 = jPart2(:)';
    jPart3 = M * ones(1, numCols);
    J = [jPart1, jPart2, jPart3];

    paddedImgZeros = zeros(numRows, numCols + P);
    paddedImgColIdxsImg = halfP + (1:numCols);
    paddedImgSIdxs = [jPart1, jPart2] - 1;
    allColsOne = ones(numCols, 1);

    j = 0;
    if isfinite(maxNumFeaturesSoughtK)
        j = maxNumFeaturesSoughtK;
    end
    foundFeatures = cell(j, 2);
    featureCount = 0;
    
    numWeightImgs = length(weightImgs);
    
    sparseMatrices = cell(length(weightImgs),1);
    
    calculateNewPath = true(2,1);
    
    for weightImgNum = 1:numWeightImgs
        paddedImg = paddedImgZeros;
        paddedImg(:, paddedImgColIdxsImg) = weightImgs{weightImgNum};
        paddedImg(:, paddedImgColIdxsNegOne) = -1;
        paddedImg = paddedImg';
        paddedImg = paddedImg(:);
        S = [paddedImg(paddedImgSIdxs); allColsOne]';

        removeIdxs = find(S == -1);
        currI = I;
        currJ = J;
        currS = S;
        currI(removeIdxs) = [];
        currJ(removeIdxs) = [];
        currS(removeIdxs) = [];

        sparseMatrices{weightImgNum} = sparse(currI, currJ, currS, M, M);
    end

    while featureCount < maxNumFeaturesSoughtK
        featurePathLens = inf(1, numWeightImgs);
        featurePaths = zeros(numRows, numWeightImgs);
        featurePathsConverted = featurePaths;
        for weightImgNum = 1:numWeightImgs
            
            if calculateNewPath(weightImgNum)
                [currFeaturePath, currFeaturePathLen, ~] = shortestpath(digraph(sparseMatrices{weightImgNum}), 1, M, 'method', 'acyclic');

%                 [currFeaturePathLen, currFeaturePath, ~] = graphshortestpath(sparseMatrices{weightImgNum}, 1, M, 'method', 'acyclic');
                if isempty(currFeaturePath)
                    featurePathLens(weightImgNum) = inf;
                    featurePaths(:, weightImgNum) = NaN;
                else
                    currFeaturePath = currFeaturePath(2:(numRows + 1));

                    featurePathLens(weightImgNum) = currFeaturePathLen;
                    featurePaths(:, weightImgNum) = currFeaturePath';

                    currFeaturePathConverted = mod(currFeaturePath - 1, numCols + P) - halfP;
                    featurePathsConverted(:, weightImgNum) = currFeaturePathConverted';
                end
            else
                featurePathLens(weightImgNum) = previousPathLen;
                featurePaths(:,weightImgNum) = previousPath;
                featurePathsConverted(:,weightImgNum) = previousPathConverted;
            end

        end


        [~,ordering] = sort(featurePathLens);
        shortPathIdx = ordering(1);
        longPathIdx = ordering(2);
        calculateNewPath(shortPathIdx) = true;
        
        if featurePathLens(shortPathIdx) == inf
            break;
        elseif isfinite(featurePathLens(longPathIdx)) && ~any(abs(featurePathsConverted(:,longPathIdx) - featurePathsConverted(:,shortPathIdx)) < trimPx+1)
            calculateNewPath(longPathIdx) = false;
            previousPathLen = featurePathLens(longPathIdx);
            previousPath = featurePaths(:,longPathIdx);
            previousPathConverted = featurePathsConverted(:,longPathIdx);
        else
            calculateNewPath(longPathIdx) = true;
        end

        featurePath = featurePaths(:,shortPathIdx);

        % Relevant pixels are now removed from both imgUp and imgDown
        for weightImgNum = 1:numWeightImgs
            sparseMatrices{weightImgNum}(1, max(featurePath(1) - trimPx, 2):min(featurePath(1) + trimPx, M - 1)) = 0;

            sparseMatrices{weightImgNum}(max(featurePath(end) - trimPx, 2):min(featurePath(end) + trimPx, M - 1), M) = 0;

            for rowNum = 1:numRows-1
                sparseMatrices{weightImgNum}(max(2, featurePath(rowNum) - trimPx - halfP):min(featurePath(rowNum) + trimPx + halfP, M - 1),...
                    max(2, featurePath(rowNum + 1) - trimPx):min(featurePath(rowNum + 1) + trimPx, M - 1)) = 0;
            end
        end

        featureCount = featureCount + 1;
        foundFeatures(featureCount, :) = {featurePathsConverted(:, shortPathIdx), featurePathLens(shortPathIdx)};  
    end

    foundFeatures = foundFeatures(1:featureCount, :);
    featurePathMeans = nan(featureCount, 1);
    for featuresNum = 1:featureCount
        featurePath = foundFeatures{featuresNum,1};
        featurePathMeans(featuresNum) = ceil(mean(featurePath));
    end
    [~, ordering] = sort(featurePathMeans);
    foundFeatures = foundFeatures(ordering, :);
end
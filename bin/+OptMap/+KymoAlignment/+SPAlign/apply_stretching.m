function [imgStretched] = apply_stretching(img, stretchFactorsMat)
    % APPLY_HORIZONTAL_STRETCHING - interpolates/stretches an image's
    %  rows with the associated horizontal stretch factors
    %
    % Inputs:
    %   img
    %     the image to stretch
    %   stretchFactorsMat
    %     matrix containing the horizontal stretch factor of each pixel
    % 
    % Outputs:
    %   imgStretched
    %     the stretched image
    % 
    % Authors:
    %   Saair Quaderi: Separated from stretch factor computation code;
    %     Refactored stretching functionality;
    %     Improved NaN handling (nearest horizontal neighbor extrapolation)
    
    imgSz = size(img);
    numRows = imgSz(1);
    numCols = imgSz(2);

    imgStretched = nan(size(img));
    matLeft = cell(1,numRows);
    for rowNum = 1:numRows
        % The aligned row is created, by linear "stretching" using the 
        % cumuative sum of the stretch factors of each pixel
        pxToStretch = find(~isnan(stretchFactorsMat(rowNum,:)));
        colIdxs = 1:length(pxToStretch);

        imgRow = img(rowNum,pxToStretch);
        matLeft{rowNum} = img(rowNum,1:pxToStretch(1)-1);
        
        oldNonnanMask = not(isnan(imgRow));
        nonnanVals = imgRow(oldNonnanMask);
        imgRowFilled = mean([...
            interp1(colIdxs(oldNonnanMask), nonnanVals, colIdxs, 'nearest', 'extrap'); ...
            fliplr(interp1(colIdxs(fliplr(oldNonnanMask)), fliplr(nonnanVals), colIdxs, 'nearest', 'extrap'))], ...
            1); % fill NaNs with the nearest nonnan value (average of nearest nonnan values if nearness is a tie)

        interpCoordsX = cumsum(stretchFactorsMat(rowNum, pxToStretch));
        stretchedImgRow = interp1(interpCoordsX, imgRowFilled, 1:max(interpCoordsX), 'linear','extrap');
        newNanMask = interp1(interpCoordsX, double(oldNonnanMask), 1:max(interpCoordsX), 'nearest') < 1;
        stretchedImgRow(newNanMask) = NaN;

        imgStretched(rowNum, 1:length(stretchedImgRow)) = stretchedImgRow;
        imgStretched(rowNum, length(stretchedImgRow)+1: length(stretchedImgRow) + length(img(rowNum,pxToStretch(end)+1:end))) = img(rowNum,pxToStretch(end)+1:end);
    end
    
    lenmat = cellfun(@(x) length(x),matLeft);
    maxLen = max(lenmat);
    leftMat = nan(numRows,maxLen);
    for i=1:numRows
        leftMat(i,end-lenmat(i)+1:end) = matLeft{i};
    end
    
    imgStretched = [leftMat imgStretched];
    
end
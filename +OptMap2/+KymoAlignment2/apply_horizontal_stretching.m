function [imgStretched] = apply_horizontal_stretching(img, stretchFactorsMat)
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

    imgStretched = img;
    colIdxs = 1:numCols;
    for rowNum = 1:numRows
        % The aligned row is created, by linear "stretching" using the 
        % cumuative sum of the stretch factors of each pixel
        imgRow = imgStretched(rowNum, :);
        oldNonnanMask = not(isnan(imgRow));
        nonnanVals = imgRow(oldNonnanMask);
        imgRowFilled = mean([...
            interp1(colIdxs(oldNonnanMask), nonnanVals, colIdxs, 'nearest', 'extrap'); ...
            fliplr(interp1(colIdxs(fliplr(oldNonnanMask)), fliplr(nonnanVals), colIdxs, 'nearest', 'extrap'))], ...
            1); % fill NaNs with the nearest nonnan value (average of nearest nonnan values if nearness is a tie)

        interpCoordsX = cumsum(stretchFactorsMat(rowNum, :));
        stretchedImgRow = interp1(interpCoordsX, imgRowFilled, 1:numCols, 'spline');
        newNanMask = interp1(interpCoordsX, double(oldNonnanMask), 1:numCols, 'nearest') < 1;
        stretchedImgRow(newNanMask) = NaN;

        imgStretched(rowNum, :) = stretchedImgRow;
    end
end
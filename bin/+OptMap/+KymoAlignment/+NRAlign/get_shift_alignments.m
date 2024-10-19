function [shiftingVect] = get_shift_alignments(unalignedKymo, maxShiftPerRow)
    % GET_SHIFT_ALIGNMENTS - Calculates the optimal shifting of each
    %   row for aligning a kymograph, without stretching it
    %
    % Inputs:
    %  unalignedKymo
    %   input image to perform shift alignment on (e.g.  a
    %    kymograph with each row representing a timeframe)
    %  maxShift
    %    the maximum shifting for each row relative to the last row
    %
    % Outputs:
    %   shiftingVect
    %    the optimal shifting of each row relative to the first row
    %    (with positive values meaning shifting right)
    %
    % Authors:
    %   Henriks Nordanger
    %   Saair Quaderi

    numRows = size(unalignedKymo, 1);
    numCols = size(unalignedKymo, 2);

    % The mean intensity is subtracted from each pixel
    unalignedKymo = unalignedKymo - mean(unalignedKymo(:));

    % Array for containing the horizontal position of each row, is
    % pre-allocated
    
    shiftingVect = zeros(1, numRows);
    shiftingVect = shiftingVect + numCols; % temporary numCols adjustment for indexing purposes

    % compareWithFirstTF = true;
    
    % Each row (except the first) is investigated
    for rowIdx = 2:numRows

        % if compareWithFirstTF
        %     comparisonRowIdx = 1;
        % else
        %     comparisonRowIdx = rowIdx - 1;
        % end
        comparisonRowIdx = 1;
        
        % The cross-correlation between row 1 and the present row is
        %  calculated
        xCorrs = xcorr(unalignedKymo(comparisonRowIdx, :), unalignedKymo(rowIdx, :));

        % The positioning of the present row, relative to previous row,
        %  is optimized to fit the row being compared to
        [~, shiftingVect(rowIdx)] = max(xCorrs(shiftingVect(rowIdx - 1) - maxShiftPerRow:shiftingVect(rowIdx - 1) + maxShiftPerRow));

        % The absolute positioning of the present row is found
        shiftingVect(rowIdx) = shiftingVect(rowIdx) + shiftingVect(rowIdx-1) - (maxShiftPerRow + 1);

    end

    shiftingVect = shiftingVect - numCols;  % remove temporary numCols adjustment for indexing purposes
end
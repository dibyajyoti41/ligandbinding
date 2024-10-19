function shiftAlignedKymo = shift_rows(unalignedKymo, shiftingVect)  
    % SHIFT_ROWS - Shifts the rows in a kymograph without stretching it.
    %
    % Inputs:
    %  unalignedKymo
    %   unaligned kymograph
    %  shiftingVect
    %   array containing the shift to be performed for each row
    %
    % Outputs:
    %  shiftAlignedKymo
    %   kymograph with rows shifted
    %
    % Authors:
    %  Henrik Nordanger
    %  Saair Quaderi
    
    padVal = NaN;
    if islogical(unalignedKymo)
        padVal = false;
    end

    numRows = size(unalignedKymo, 1);
    numCols = size(unalignedKymo, 2);
    doubleNumCols = numCols * 2;

    %Check if the two input arguments are compatible
    
    if length(shiftingVect) ~= numRows
        disp('The number of specified row shifts do not match the number of rows');
    elseif (max(shiftingVect) > numCols) || (min(shiftingVect) < -numCols)
        disp('Row shift values are outside acceptable range');
    else
        % NaN values are added at the sides of the kymograph, to make
        %  room for any shifted row
        unalignedKymo = padarray(unalignedKymo, [0, numCols], padVal);

        % Each row (except the first) is investigated.
        for rowNum = 2:numRows
            % The present row is temporarily and separately saved
            tempRow = unalignedKymo(rowNum, (numCols + 1):doubleNumCols);

            % The whole row is set to NaN
            unalignedKymo(rowNum, numCols+1:doubleNumCols) = repmat(padVal, [1, numCols]);

            % The temporarily saved row is placed where specified by
            %  therelevant element of 'shifts'
            shift = shiftingVect(rowNum);
            unalignedKymo(rowNum, (numCols + 1 + shift):(doubleNumCols + shift)) = tempRow;
        end

        %The relevant part of the resulting kymograph is returned
        shiftAlignedKymo = unalignedKymo(:,(min(shiftingVect) + numCols + 1):(max(shiftingVect) + doubleNumCols));
    end
end
function barcode = normalize_intensity(barcode, barcodeBitmask, setMinAsZeroMaxAsOne)
    % Shifts and scales the intensities for matrix A between 0 and 1
    %  If setMinAsZeroMaxAsOne is set to true:
    %   all values are shifted and scaled linearly such that
    %   the minimum value is moved to 0
    %   the maximum value is moved to 1
    %   (Unless all values are equal in which case everything is set to 0.5)
    %  Otherwise:
    %   an intensity of 0 is always shifted to 0.5 and
    %   the maximum value is at most 1 and minimum value is at least 0
    %   and either the maximum is 1, the minimum is 0, or all are 0.5
    if nargin < 2
        barcodeBitmask = true(size(barcode));
    end

    if nargin < 3
        setMinAsZeroMaxAsOne = false;
    end

    minVal = min(barcode(barcodeBitmask));
    maxVal = max(barcode(barcodeBitmask));
    if setMinAsZeroMaxAsOne
        valRange = maxVal - minVal;
        if (valRange > 0)
            barcode = (barcode - minVal); % shift such that minimum is zero
            barcode = barcode / valRange; % scale such that maximum is one
        else
            barcode(:) = 1 / 2;
        end
    else % shifting so 0 is 0.5 and scaling to fit interval
        scaling = max(abs(minVal), abs(maxVal));
        if (scaling > 0) % if not all values are 0
            % linearly scale (without shifting zero)
            %  such that all values are within [-1, 1]
            %  but as spread apart as possible
            barcode = barcode / scaling;
        end
        % shift and linearly scale values from [-1, 1] to [0, 1]
        barcode = (1 + barcode) / 2;
    end
    barcode = max(0, min(1, barcode));
end
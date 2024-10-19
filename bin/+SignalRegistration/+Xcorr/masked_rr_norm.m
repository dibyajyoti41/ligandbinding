function [seq] = masked_rr_norm(seq, bitmask, zeroifyRest)
    % MASKED_RR_NORM - Scales a barcode such that the portion of
    %  the curve corresponding to bitmask values of 1 have a
    %  mean = 0 and std dev = 1
    %  (bitmasked version of reisner-rescaling a.k.a. zscore)
    %  
    % Inputs:
    %  seq (vector of the data sequence)
    %  bitmask (the bitmask with 1 at the indices of the barcode
    %     where the value should be taken into account for rescaling)
    %  zeroifyRest (optional boolean, which, if true, will set all the
    %     rescaled barcode values at indices where the bitmask
    %     is 0 to a value of 0)
    %
    % Outputs:
    %   seq (the rescaled sequence)
    % 
    % Dependencies: built-in matlab functions
    %
    % By: Saair Quaderi

    if isempty(seq) && isempty(bitmask)
        return;
    end
    validateattributes(seq, {'numeric'}, {'vector'});
    validateattributes(bitmask, {'logical', 'numeric'}, {'vector', 'binary'});
    if (length(bitmask) ~= length(seq))
        error('Bad bitmask length for reisner rescaling normalization');
    end
    if nargin < 3
        zeroifyRest = false;
    else
        validateattributes(zeroifyRest, {'logical', 'numeric'}, {'scalar', 'binary'});
    end
    if not(islogical(bitmask))
        bitmask = logical(bitmask);
    end
    subset = seq(bitmask);
    if not(isempty(subset))
        seq = (seq - mean(subset));
        variance = std(subset);
        if (variance ~= 0)
            seq = seq/variance;
        end
    end
    if zeroifyRest
        seq(~bitmask) = 0;
    end
end

function [extendedSeq] = extend_circular_seq(seqPeriod, lenPaddingFront, lenPaddingEnd)
    % EXTEND_CIRCULAR_SEQ - Extends a circular sequence by padding
    %  it by a certain number of values before/after it
    %  (the padded values are determined by the seq)
    %  
    % Inputs:
    %  seqPeriod (period of the sequence)
    %  lenPaddingFront (the amount to extend the sequence at the
    %    front)
    %  lenPaddingEnd (the amount to extend the sequence at the
    %    front)
    %
    % Outputs:
    %  extendedSeq (the extended sequence)
    % 
    % Dependencies: built-in matlab functions
    %
    % By: Saair Quaderi

    validateattributes(lenPaddingFront, {'numeric'}, {'scalar', 'integer', 'nonnegative'});
    validateattributes(lenPaddingEnd, {'numeric'}, {'scalar', 'integer', 'nonnegative'});
    validateattributes(seqPeriod, {'numeric', 'logical'}, {'column', 'nonempty'});
    if (lenPaddingFront == 0) && (lenPaddingEnd == 0)
        extendedSeq = seqPeriod;
        return;
    end
    asLogical = islogical(seqPeriod);
    lenSeq = length(seqPeriod);
    frontPadding = repmat(seqPeriod, ceil(lenPaddingFront/lenSeq), 1);
    endPadding = repmat(seqPeriod, ceil(lenPaddingEnd/lenSeq), 1);
    extendedSeq = [frontPadding((end + 1 - lenPaddingFront):end); seqPeriod; endPadding(1:lenPaddingEnd)];
    if (asLogical)
        extendedSeq = logical(extendedSeq);
    end
end
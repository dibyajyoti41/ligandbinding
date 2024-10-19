function [extendedSeq] = zero_extend_linear_seq(seqPeriod, lenPaddingFront, lenPaddingEnd)
    % ZERO_EXTEND_LINEAR_SEQ - Extends a linear sequence by padding
    %  it by a certain number of zeros before/after it
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
    if not(isempty(seqPeriod))
        validateattributes(seqPeriod, {'numeric', 'logical'},{'column'});
    else
        validateattributes(seqPeriod, {'numeric', 'logical'}, {});
    end
    if (lenPaddingFront == 0) && (lenPaddingEnd == 0)
        extendedSeq = seqPeriod;
        return;
    end
    asLogical = islogical(seqPeriod);
    extendedSeq = [zeros(lenPaddingFront, 1); seqPeriod; zeros(lenPaddingEnd, 1)];
    if (asLogical)
        extendedSeq = logical(extendedSeq);
    end
end

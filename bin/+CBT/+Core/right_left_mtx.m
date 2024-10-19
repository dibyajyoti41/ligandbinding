function [ tmpNormConstants, tmpConstraintMat] = right_left_mtx_ml( seqStart, interval, seqEnd, tmpBindingConstants, tmpBindingIdxs, tmpVectorPrev, tmpTransferMatrix, mtimesReverse)
%RIGHT_LEFT_MTX Summary of this function goes here
%   Detailed explanation goes here
     tmpNormConstants(seqEnd) = 0;
     tmpConstraintMat(9, seqEnd) = 0;
    for seqBpIdx = seqStart:interval:seqEnd
        tmpTransferMatrix(tmpBindingIdxs) = tmpBindingConstants(:, seqBpIdx);
        if mtimesReverse
            tmpVectorCurr = mtimes(tmpTransferMatrix, tmpVectorPrev);
        else
            tmpVectorCurr = mtimes(tmpVectorPrev, tmpTransferMatrix);
        end
            
        tmpValueCurrNorm = norm(tmpVectorCurr);
        tmpVectorCurr = tmpVectorCurr / tmpValueCurrNorm;

        tmpNormConstants(seqBpIdx) = tmpValueCurrNorm;
        tmpConstraintMat(:, seqBpIdx) = tmpVectorCurr;
        tmpVectorPrev = tmpVectorCurr;
    end
end


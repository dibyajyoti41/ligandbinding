function [ bindingProb ] = binding_probs_mtx_ml( numBasepairs, seqBindingConstants, bindingIdx, transferMatrix, firstVec, lastVec, rightMatrix, leftMatrix )
%BINDING_PROBS_MTX_ML Summary of this function goes here
%   Detailed explanation goes here

        bindingProb = zeros(1, numBasepairs);
        seqBpIdx = 1;
        transferMatrix(bindingIdx) = seqBindingConstants(seqBpIdx);
        bindingProb(seqBpIdx) = mtimes(mtimes(firstVec, transferMatrix), rightMatrix(:, seqBpIdx + 1));
        for seqBpIdx = 2:(numBasepairs - 1)
            transferMatrix(bindingIdx) = seqBindingConstants(seqBpIdx);
            bindingProb(seqBpIdx) = mtimes(mtimes(leftMatrix(:, seqBpIdx - 1)', transferMatrix), rightMatrix(:, seqBpIdx + 1));
        end
        seqBpIdx = numBasepairs;
        transferMatrix(bindingIdx) = seqBindingConstants(seqBpIdx);
        bindingProb(seqBpIdx) = mtimes(mtimes(leftMatrix(:, seqBpIdx - 1)', transferMatrix), lastVec);

end


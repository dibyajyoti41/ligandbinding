function [seqBindingConstants] = get_sequence_binding_constants(bindingConstantsACGT, ntBitsmart_theorySeq_uint8, showWarning, defaultBindingConstant)
    % GET_SEQUENCE_BINDING_CONSTANTS - Gets the local sequence-dependent
    %  binding constants associated with each basepair along the nucleotide
    %  sequence
    %
    % Inputs:
    %  bindingConstantsACGT
    %    4x4x...x4 matrix with k dimensions where k is the number of bps
    %    occupied by the ligand and each value is associated with each of
    %    the possible k-basepair sequences (ordering in each dimension is
    %    A, C, G, T and jth dimension represents the jth basepair in the
    %    sequence)
    %  ntBitsmart_theorySeq_uint8
    %    dna nucleotide label sequence in our custom uint8 bitsmart
    %    encoding
    %  showWarning (optional, defaults to true)
    %    shows a warning if there were ambiguous basepair labels that yield
    %    multiple possible binding constant results which had to be
    %    averaged
    %  defaultBindingConstant (optional, defaults to 0)
    %    filler binding constant value for end of sequence since
    %    there can be no overhang and circularity is not implied at this
    %    level of the code. In theory, this shouldn't actually make a
    %    difference to results.
    %
    % Outputs:
    %  seqBindingConstants
    %    the binding constants associated with each of the basepairs
    %
    % Authors:
    %   Saair Quaderi
    
    import NtSeq.Core.uint8_in_binary;

    if nargin < 3
        showWarning = true;
    end

    if nargin < 4
        defaultBindingConstant = 0;
    end
    bpsOccupied = ndims(bindingConstantsACGT);
    seqLen = length(ntBitsmart_theorySeq_uint8);
    bpSeqIdxHelper = floor(log2(double(ntBitsmart_theorySeq_uint8)));
    ambiguousIdxs = find(uint8(2.^bpSeqIdxHelper) ~= ntBitsmart_theorySeq_uint8);
    ambiguousIdxsExist = not(isempty(ambiguousIdxs));
    ambiguousIdxsExtendedBack = [];
    if ambiguousIdxsExist
        % extend indexes back to include all start indices which
        %  will yield results from an ambiuous labeling
        ambiguousIdxsExtendedBack = bsxfun(@plus, (1 - bpsOccupied):0, ambiguousIdxs(:));
        ambiguousIdxsExtendedBack = unique(ambiguousIdxsExtendedBack(:));
        ambiguousIdxsExtendedBack = ambiguousIdxsExtendedBack(ambiguousIdxsExtendedBack > 0);

        % extend forward to include all indices looked at from any
        %  start index that will include any ambiguous labeling
        ambiguousIdxsExtendedFull = bsxfun(@plus, 0:(bpsOccupied - 1), ambiguousIdxsExtendedBack(:));
        ambiguousIdxsExtendedFull = unique(ambiguousIdxsExtendedFull(:));
        ambiguousIdxsExtendedFull = ambiguousIdxsExtendedFull(ambiguousIdxsExtendedFull <= seqLen);
    end

    seqBindingConstantsIdxs = round(conv(bpSeqIdxHelper, 4.^(0:(bpsOccupied - 1)), 'valid')) + 1;
    seqBindingConstants = bindingConstantsACGT(seqBindingConstantsIdxs);
    seqBindingConstants = [seqBindingConstants(:)', repmat(defaultBindingConstant, 1, bpsOccupied - 1)];
    if ambiguousIdxsExist
        seqBindingConstants(ambiguousIdxsExtendedBack) = NaN;

        bpSeq_logical = false(seqLen, 4);
        bpSeq_logical(ambiguousIdxsExtendedFull, :) = uint8_in_binary(ntBitsmart_theorySeq_uint8(ambiguousIdxsExtendedFull), 4);

        runsAverageCount = 0;
        for idxIdx = 1:numel(ambiguousIdxsExtendedBack)
            seqIdxBindingStart = ambiguousIdxsExtendedBack(idxIdx);
            if (seqIdxBindingStart + bpsOccupied) > seqLen
                continue;
            end
            bpSubseqLogicalIdxs = bpSeq_logical(seqIdxBindingStart - 1 + (1:bpsOccupied), :);
            bpSubseqLogicalIdxs = mat2cell(bpSubseqLogicalIdxs(1:bpsOccupied, :), ones(1, bpsOccupied), 4);
            bpBindingConstant = bindingConstantsACGT(bpSubseqLogicalIdxs{:});
            bpBindingConstant = bpBindingConstant(:);
            if any(bpBindingConstant ~= bpBindingConstant(1))
                runsAverageCount = runsAverageCount + 1;
                seqBindingConstants(seqIdxBindingStart) = mean(bpBindingConstant);
            else
                seqBindingConstants(seqIdxBindingStart) = bpBindingConstant(1);
            end
        end
        if showWarning && (runsAverageCount > 0)
            warning('Multiple different binding constants associated with ambiguous nucleotide labeling in the input sequence are being averaged');
        end
    end
end

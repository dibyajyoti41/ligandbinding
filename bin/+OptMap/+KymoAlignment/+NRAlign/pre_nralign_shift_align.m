function [shiftAlignedKymo, shiftingVect] = pre_nralign_shift_align(unalignedKymo, rowSmoothingWindowLen_pixels, blurSigmaWidth_pixels)
    % PRE_NRALIGN_SHIFT_ALIGN - non-recursively aligns a kymograph,
    %  using a modified wpalign algorithm.
    %
    % Inputs:
    %	unalignedKymo
    %	  the unaligned kymograph
    %
    % Outputs:
    %   shiftAlignedKymo
    %     the shift aligned kymograph on which stretch factors were
    %     computed to produce the aligned kymograph
    %   shiftingVect
    %     the amount each row was shifted
    %
    % Authors:
    %	Henrik Nordanger
    %   Saair Quaderi


    % Each row of the kymograph is smoothed out in order to reduce noise
    % SQ: Is this necessary? are these parameters good choices?
    if nargin < 2
        rowSmoothingWindowLen_pixels = 10;
    end
    if nargin < 3
        blurSigmaWidth_pixels = 2;
    end

    rowBlurredUnalignedKymo = unalignedKymo;
    if rowSmoothingWindowLen_pixels ~= 1
        hSize = [1, rowSmoothingWindowLen_pixels];
        blurringKernel = fspecial('gaussian', hSize, blurSigmaWidth_pixels);
        rowBlurredUnalignedKymo = conv2(rowBlurredUnalignedKymo, blurringKernel, 'valid');
    end

    %The kymograph is "prealigned" using shift-based alignment
    % The maximum shifting of each row
    maxShiftPerRow = 3;
    import OptMap.KymoAlignment.NRAlign.get_shift_alignments;
    shiftingVect = get_shift_alignments(rowBlurredUnalignedKymo, maxShiftPerRow);

    import OptMap.KymoAlignment.NRAlign.shift_rows;
    shiftAlignedKymo = shift_rows(unalignedKymo, shiftingVect);
end
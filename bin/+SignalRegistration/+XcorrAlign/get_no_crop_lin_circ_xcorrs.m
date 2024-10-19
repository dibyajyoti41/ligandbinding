function [xcorrs, coverageLens, firstOffset] = get_no_crop_lin_circ_xcorrs(linSeq, circSeq, linSeqBitmask, circSeqBitmask)
    % get_no_crop_lin_circ_xcorrs
    %  Computes the normalized cross-correlation for a linear
    %   sequence and a circular sequence, flipped/not flipped,
    %   at different circular shifts (but without any cropping of
    %   the linear sequence provided)
    % Inputs:
    %  linSeq & circSeq - the non-periodic and periodic sequence
    %    we wish to align
    %  linSeqBitmask & circSeqBitmask - bitmasks for the
    %   two sequences
    %   (linSeqBitmask & circSeqBitmask are aligned based only on
    %     the values and positions of the true bits in the bitmask
    %     corresponding to the data in the sequences)
    %
    % Outputs:
    %  xcorrs - raw Pearson cross-correlation coefficient
    %    values as computed in a three-dimensional array.
    %    The first dimensions specifies whether a sequence was
    %     flipped (index 2 if true, 1 if false).
    %    The second dimension specifies how much circular
    %      shifting there was (index 1 = no shift, 2 = shift of 1,
    %    etc.)
    %    The third dimension typically specifies the linear offset,
    %      (index k => firstOffset + k linear shifting)
    %      but is not used for any offsets in this case
    %  coverageLens - Effective sample sizes corresponding to each
    %    xcorrs value
    import SignalRegistration.XcorrAlign.norm_xcorr;

    linSeq = linSeq(:)';
    circSeq = circSeq(:)';
    if nargin < 3
        linSeqBitmask = true(size(linSeq));
    end
    if nargin < 4
        circSeqBitmask = true(size(circSeq));
    end
    linSeqBitmask = linSeqBitmask(:)';
    circSeqBitmask = circSeqBitmask(:)';
    firstOffset = 0;
    circSeqLen = length(circSeq);
    linSeqLen = length(linSeq);
    lenDiff = circSeqLen - linSeqLen;
    if lenDiff < 0
        error('Linear sequence was longer than circular sequence.');
    end
    zPadding = zeros(1, lenDiff);
    flippedLinSeqPadded = [zPadding, fliplr(linSeq)];
    flippedLinSeqPaddedBitmask = [logical(zPadding), fliplr(linSeqBitmask)];

    [xcorrsReg, coverageLensReg, ~] = norm_xcorr(linSeq, circSeq, linSeqBitmask, circSeqBitmask, false, true);
    [xcorrsFlipped, coverageLensFlipped, ~] = norm_xcorr(flippedLinSeqPadded, circSeq, flippedLinSeqPaddedBitmask, circSeqBitmask, false, true);
    xcorrs = [
        xcorrsReg;...
        xcorrsFlipped...
    ];
    coverageLens = [
        coverageLensReg;...
        coverageLensFlipped...
    ];
    colIndices = (((size(xcorrs, 2) + 1)/2) - (0:(circSeqLen - 1)));
    xcorrs = xcorrs(:, colIndices);
    coverageLens = coverageLens(:, colIndices);
end
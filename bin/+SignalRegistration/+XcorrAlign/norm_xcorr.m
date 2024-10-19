function [xcorrs, coverageLens, delays] = norm_xcorr(a, b, bitmaskA, bitmaskB, circularA, circularB, delay)
    % norm_xcorr
    %  Quickly compute Pearson cross-correlation coefficients
    %   of two discrete bitmasked sequences after normalizing their
    %   effective samples to a mean of 0 and variance of 1
    %   (where the effective samples are determined by ANDing their
    %    bitmasks at each delay)
    %
    % Inputs:
    %  seqA & seqB - the two sequences we wish to compute
    %   pearson cross-correlations coefficients for
    %  bitmaskA & bitmaskB - bitmasks for seqA and seqB
    %   (seqA & seqB are aligned based only on the values
    %     and positions of the true bits in the bitmask
    %     corresponding to the data in the sequences_
    %  circularA & circularB - whether the sequences
    %   given should be treated like they are circular
    %   (e.g. allow circularly shifting them, respecting their
    %     periodic nature) as opposed to linear sequences
    %    (which cannot be circularly shifted and don't repeat)
    %  delay - the magnitude of delay (shifting/offsetting)
    %    permitted
    %
    %  xcorrs - Raw Pearson cross-correlation coefficient
    %    values as a vector with values at different delays
    %  coverageLens - Effective sample sizes corresponding to each
    %    xcorrs value
    %  delays - The amount of delay associated with each value
    %    in xcorrs and coverageLens
    import SignalRegistration.Xcorr.masked_norm_xcorr;

    a = zscore(a(:)); % normalize
    b = zscore(b(:)); % normalize

    if nargin < 3
        bitmaskA = true(size(a));
    else
        bitmaskA = bitmaskA(:);
    end
    if nargin < 4
        bitmaskB = true(size(b));
    else
        bitmaskB = bitmaskB(:);
    end
    if nargin < 5
        circularA = false;
    end
    if nargin < 6
        circularB = false;
    end

    maxLen = max(length(a), length(b));
    if nargin < 7
        delay = max(0, maxLen - 1);
    else
        delay = abs(delay);
    end
    delay = [-1*delay, delay];

    [xcorrs, delays, coverageLens] = masked_norm_xcorr(a, b, bitmaskA, bitmaskB, circularA, circularB, delay, false);
    xcorrs = xcorrs';
    coverageLens = coverageLens';
end
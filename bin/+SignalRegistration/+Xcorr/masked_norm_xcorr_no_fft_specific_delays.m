function [xcorrs, specificDelays, coverageLens] = masked_norm_xcorr_no_fft_specific_delays(a, b, aBitmask, bBitmask, bIsCircular, specificDelays)
    % MASKED_NORM_XCORR_NO_FFT_SPECIFIC_DELAYS - Computes pearson
    %   cross-correlations coefficients for specific delays within
    %   a delay range without using fourier transforms
    %  
    % Inputs:
    %  a & b (column vectors representing the
    %    sequences to compute the cross-correlation over --
    %    note: b might represent periodic data)
    %  aBitmask & bBitmask (column vector bitmasks representing the
    %    regions on a and b which should be counted towards the-
    %    sliding dot products with 1s)
    %  bIsCircular (true if b represents a period from a repeating
    %    sequence -- will cause it to be treated as if it is
    %    indefinitely preceded and followed by its own sequence,
    %    otherwise it will be treated as if it is infinitely
    %    preceded and followed by its mean value)
    %  specificDelays (vector of non-positive values representing
    %     the delays which are of interest)
    %
    % Outputs:
    %  xcorrs (column vector containing normalized
    %    cross-correlations -- a.k.a. normalized sliding dot
    %    products -- for the specified delays)
    %  specificDelays (row vector containing the delay which
    %    corresponds to each index in xcorrs)
    import SignalRegistration.Xcorr.extend_circular_seq;
    import SignalRegistration.Xcorr.zero_extend_linear_seq;
    import SignalRegistration.Xcorr.masked_rr_norm;

    validateattributes(a, {'numeric', 'logical'}, {'column', 'nonempty', 'real', 'finite', 'nonnan'});
    validateattributes(b, {'numeric', 'logical'}, {'column', 'nonempty', 'real', 'finite', 'nonnan'});
    validateattributes(aBitmask, {'logical', 'numeric'}, {'vector', 'binary'});
    validateattributes(bBitmask, {'logical', 'numeric'}, {'vector', 'binary'});
    validateattributes(bIsCircular, {'logical', 'numeric'}, {'scalar', 'binary'});

    if isempty(specificDelays)
        xcorrs = zeros(0, 1);
        return;
    end
    validateattributes(specificDelays, {'numeric'}, {'integer', 'vector'});

    numDelays = length(specificDelays);
    xcorrs = zeros(numDelays, 1);
    coverageLens = zeros(numDelays, 1);
    minDelay = min(0, min(specificDelays));
    absMinDelay = abs(minDelay);
    maxDelay = max(0, max(specificDelays));

    n = length(b);
    bDelayZeroIndex = maxDelay + 1;
    if (bIsCircular)
        bExtended = extend_circular_seq(b, bDelayZeroIndex - 1, absMinDelay);
        bBitmaskExtended = extend_circular_seq(bBitmask, bDelayZeroIndex - 1, absMinDelay);
    else
        bExtended = zero_extend_linear_seq(b, bDelayZeroIndex - 1, absMinDelay);
        bBitmaskExtended = zero_extend_linear_seq(bBitmask, bDelayZeroIndex - 1, absMinDelay);
    end
    delayBExtendedIndicesPlusDelay = (1:n) - 1 + bDelayZeroIndex;
    for delayIdx = 1:numDelays
        delayBIndices = delayBExtendedIndicesPlusDelay - specificDelays(delayIdx);
        delayBitmask = aBitmask & bBitmaskExtended(delayBIndices);
        delaySampleLen = sum(delayBitmask);
        if (delaySampleLen > 1)
            aSampleNormalized = masked_rr_norm(a, delayBitmask, true);
            bDelaySampleNormalized = masked_rr_norm(bExtended(delayBIndices), delayBitmask, true);

            delayXcorr = sum(aSampleNormalized.*bDelaySampleNormalized);
            if (delayXcorr ~= 0) % Avoid 0/0=NaN, keep as 0
                delayXcorr = delayXcorr/(delaySampleLen - 1);
            end
        else
            delayXcorr = 0;
        end
        coverageLens(delayIdx) = delaySampleLen;
        xcorrs(delayIdx) = delayXcorr;
    end
end

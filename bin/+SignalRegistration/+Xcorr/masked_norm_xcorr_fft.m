function [xcorrs, delays, coverageLens] = masked_norm_xcorr_fft(a, b, aBitmask, bBitmask, bIsCircular, minDelay, maxDelay)
    % MASKED_NORM_XCORR_FFT - Computes pearson cross-correlations
    %  coefficients within a delay range using fourier transforms
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
    %  minDelay (non-positive value representing the minimum delay
    %    we are interested in for the range of delays)
    %  maxDelay (non-negative value representing the maximum delay
    %    we are interested in for the range of delays)
    %
    % Outputs:
    %  xcorrs (column vector containing normalized
    %    cross-correlations  -- a.k.a. normalized sliding dot
    %    products -- for the delay range)
    %  delays (row vector containing the delay which
    %    corresponds to each index in sdps)
    import SignalRegistration.Xcorr.circ_sdps_fft;
    import SignalRegistration.Xcorr.lin_sdps_fft;
    import SignalRegistration.Xcorr.masked_norm_xcorr_no_fft_specific_delays;

    validateattributes(a, {'numeric', 'logical'}, {'column', 'nonempty', 'real', 'finite', 'nonnan'});
    validateattributes(b, {'numeric', 'logical'}, {'column', 'nonempty', 'real', 'finite', 'nonnan'});
    validateattributes(aBitmask, {'logical', 'numeric'}, {'vector', 'binary'});
    validateattributes(bBitmask, {'logical', 'numeric'}, {'vector', 'binary'});
    validateattributes(bIsCircular, {'logical', 'numeric'}, {'scalar', 'binary'});

    lenA = length(aBitmask);
    lenB = length(bBitmask);
    n = max(lenA, lenB);

    if (length(aBitmask) ~= lenA)
        error('Bitmask length does not match sequence length of a')
    end
    if (length(bBitmask) ~= lenB)
        error('Bitmask length does not match sequence length of b')
    end

    if nargin < 6
        minDelay = 1 - n;
    else
        validateattributes(minDelay, {'numeric'}, {'scalar', 'integer', '<=', 0});
    end
    if nargin < 7
        maxDelay = n - 1;
    else
        validateattributes(maxDelay, {'numeric'}, {'scalar', 'integer', 'nonnegative'});
    end

    if bIsCircular
        xcorr_fft = @(refSeq, circSeq) circ_sdps_fft(refSeq, circSeq, minDelay, maxDelay);
    else
        xcorr_fft = @(refSeq, circSeq) lin_sdps_fft(refSeq, circSeq, minDelay, maxDelay);
    end

    % Normalize values initially to avoid any risk of "catastrophic cancellation"
    %  issues (todo: investigate - is this really necessary?)
    a = zscore(a);
    b = zscore(b);

    delays = minDelay:maxDelay;
    a_aBitmask = aBitmask.*a;
    b_bBitmask = bBitmask.*b;
    aa_aBitmask = a_aBitmask.*a;
    bb_bBitmask = b_bBitmask.*b;

    delaySampleSizes = xcorr_fft(aBitmask, bBitmask);
    delaySampleSizes = round(delaySampleSizes); %correct rounding errors -- make integers (should be perfectly precise)

    delaySampleASums = xcorr_fft(a_aBitmask, bBitmask);
    delaySampleAMeans = delaySampleASums./delaySampleSizes;

    delaySampleBSums = xcorr_fft(aBitmask, b_bBitmask);
    delaySampleBMeans = delaySampleBSums./delaySampleSizes;

    delayXcorrNumerators = xcorr_fft(a_aBitmask, b_bBitmask);
    delayXcorrSubtractFromNumerators = delaySampleSizes.*delaySampleAMeans.*delaySampleBMeans;
    delayXcorrNumerators = delayXcorrNumerators - delayXcorrSubtractFromNumerators;

    delayXcorrDenominatorRootFactorA = xcorr_fft(aa_aBitmask, bBitmask);
    delayXcorrSubtractFromDenominatorRootFactorA = delaySampleSizes.*(delaySampleAMeans.^2);
    delayXcorrDenominatorRootFactorA = delayXcorrDenominatorRootFactorA - delayXcorrSubtractFromDenominatorRootFactorA;

    delayXcorrDenominatorRootFactorB = xcorr_fft(aBitmask, bb_bBitmask);
    delayXcorrSubtractFromDenominatorRootFactorB = delaySampleSizes.*(delaySampleBMeans.^2);
    delayXcorrDenominatorRootFactorB = delayXcorrDenominatorRootFactorB - delayXcorrSubtractFromDenominatorRootFactorB;

    delayXcorrDenominatorRoots = max(delayXcorrDenominatorRootFactorA.*delayXcorrDenominatorRootFactorB, 0); % otherwise tiny negatives from rounding errors can make square roots imaginary
    delayXcorrDenominators = sqrt(delayXcorrDenominatorRoots);

    % When denominators are lower than 1, cross-correlations are
    % amplified, but because they are so small they can be heavily
    % affected by approximation/rounding errors from fft
    % sliding-dot-product calculations and floating point math
    % therefore we want to ignore these values and use slower more
    % precise means to get their values later.
    % In practice with normalized data this only occurs when
    % the effective sample size at a delay is very small
    supplementIndices = (delayXcorrDenominators < 1);

    % we know they should have xcorr of 0 and since they tend to
    % have bad rounding errors or divide by zero issues
    % we can just deal with them independently to avoid
    % NaN, -Inf, Inf
    makeZeroIndices = (delaySampleSizes < 2)|(delayXcorrDenominators == 0)|(delayXcorrNumerators == 0); 

    delayXcorrNumerators(makeZeroIndices) = 0;
    delayXcorrDenominators(makeZeroIndices) = 1;

    xcorrs = delayXcorrNumerators./delayXcorrDenominators;

    % supplement values as explained earlier
    xcorrs_supplement = masked_norm_xcorr_no_fft_specific_delays(a, b, aBitmask, bBitmask, bIsCircular, delays(supplementIndices));
    xcorrs(supplementIndices) = xcorrs_supplement;

    xcorrs = min(1, max(-1, xcorrs)); % correct rounding errors that lead to values exceeding -1:1 range
    coverageLens = delaySampleSizes;
end

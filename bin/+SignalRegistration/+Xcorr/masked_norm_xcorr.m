function [xcorrs, delays, coverageLens] = masked_norm_xcorr(a, b, aBitmask, bBitmask, aIsCircular, bIsCircular, delay, noFFT)
    % MASKED_NORM_XCORR - Computes pearson cross-correlations
    %   coefficients within a delay range
    %  
    % Inputs:
    %  a & b (column or row vectors representing the
    %    sequences to compute the cross-correlation over --
    %    note: a & b may be circular but if they are both
    %    circular they must be of the same length and if
    %    only one is circular it must not have shorter length
    %    than the other one)
    %  aBitmask & bBitmask (column vector bitmasks representing the
    %    regions on a and b which should be counted towards the-
    %    sliding dot products with 1s)
    %  aIsCircular & bIsCircular (true respectively when a & b
    %    represent a period from a repeating sequence -- will cause
    %    it to be treated as if it is indefinitely preceded and
    %    followed by its own sequence, otherwise each will be
    %    treated as if it is infinitely preceded and followed by
    %    its mean value)
    %  delay (delay range [minDelay, maxDelay] or scalar value
    %     such that delay range is [-abs(delay), abs(delay)]
    %     note: minDelay greater than 0 is not supported
    %     note: maxDelay less than 0 is not supported)
    %  noFFT (don't use fast fourier transforms to calculate value)
    %
    % Outputs:
    %  xcorrs (vector containing normalized pearson
    %    cross-correlations -- a.k.a. normalized sliding dot
    %    products -- for the specified delays)
    %  delays (row vector containing the delays which
    %    corresponds to each index in xcorrs)
    import SignalRegistration.Xcorr.masked_norm_xcorr_no_fft_specific_delays;
    import SignalRegistration.Xcorr.zero_extend_linear_seq;
    import SignalRegistration.Xcorr.unvalidated_masked_norm_xcorr;

    if isempty(a) && not((nargin >= 3) && not(isempty(aBitmask)))
        a = 0;
        aBitmask = false;
    end

    if isempty(b) && not((nargin >= 4) && not(isempty(bBitmask)))
        b = 0;
        bBitmask = false;
    end

    validateattributes(a, {'numeric'}, {'nonempty', 'vector', 'real', 'finite', 'nonnan'}, 1);
    validateattributes(b, {'numeric'}, {'nonempty', 'vector', 'real', 'finite', 'nonnan'}, 2);

    if nargin < 3
        aBitmask = true(size(a));
    else
        validateattributes(aBitmask, {'logical', 'numeric'}, {'vector', 'binary', 'size', size(a)}, 3);
        if not(islogical(aBitmask))
            aBitmask = logical(aBitmask);
        end
    end
    if nargin < 4
        bBitmask = true(size(b));
    else
        validateattributes(bBitmask, {'logical', 'numeric'}, {'vector', 'binary', 'size', size(b)}, 4);
        if not(islogical(bBitmask))
            bBitmask = logical(bBitmask);
        end
    end

    asRow = isrow(a) && isrow(b);

    % Represent as column vectors and normalize to prevent risk
    %  of catastrophic cancelations
    a = a(:);
    b = b(:);
    aBitmask = aBitmask(:);
    bBitmask = bBitmask(:);
    lenA = length(a);
    lenB = length(b);

    if nargin < 5
        aIsCircular = false;
    else
        validateattributes(aIsCircular, {'logical', 'numeric'}, {'scalar', 'binary'}, 5);
        aIsCircular = logical(aIsCircular);
    end

    if nargin < 6
        bIsCircular = false;
    else
        validateattributes(bIsCircular, {'logical', 'numeric'}, {'scalar', 'binary'}, 6);
        bIsCircular = logical(bIsCircular);
    end

    maxLen = max(lenA, lenB);
    % minLag must be some finite non-positive integer
    % maxLag must be some finite non-negative integer
    if (nargin < 7) || isempty(delay)
        delay = maxLen - 1;
        minDelay = -1*delay;
        maxDelay = delay;
    else
        if not(isscalar(delay))
            validateattributes(delay, {'numeric'}, {'integer', 'vector', 'size', [1,2], 'nondecreasing'}, 7);
            minDelay = delay(1);
            maxDelay = delay(2);
            if minDelay > 0
                error('The minDelay, delay(1), must be non-positive');
            end
            if maxDelay < 0
                error('The maxDelay, delay(2), must be non-negative');
            end
        else
            validateattributes(delay, {'numeric'}, {'integer'}, 7);
            delay = abs(delay);
            minDelay = -1*delay;
            maxDelay = delay;
        end
    end

    if nargin < 8
        noFFT = false;
    else
        validateattributes(noFFT, {'logical', 'numeric'}, {'scalar', 'binary'}, 8);
        noFFT = logical(noFFT);
    end

    if (aIsCircular && bIsCircular && (lenA ~= lenB))
        error('Unable to compute cross-correlation coefficients for circular sequences of differing lengths');
    elseif ((bIsCircular && (lenB < lenA)) || (aIsCircular && (lenA < lenB)))
        error('Unable to compute cross-correlation coefficients when circular sequence is shorter than linear sequence');
    elseif (aIsCircular && (lenB < lenA))
        % swap a and b
        % the base functionality assumes that either both are
        % linear sequences or b is a circular sequence
        % and is not shorter than a
        [xcorrs, delays] = masked_norm_xcorr_no_fft_specific_delays(b, a, bBitmask, aBitmask, bIsCircular, aIsCircular, [-maxDelay, -minDelay], noFFT);
        % correct the lags since they were swapped
        delays = -fliplr(delays);
        return;
    elseif not(aIsCircular || bIsCircular) && (lenA ~= lenB)
        if lenA < maxLen
            a = zero_extend_linear_seq(a, 0, maxLen - lenA);
            aBitmask = zero_extend_linear_seq(aBitmask, 0, maxLen - lenA);
        else
            b = zero_extend_linear_seq(a, 0, maxLen - lenB);
            bBitmask = zero_extend_linear_seq(bBitmask, 0, maxLen - lenB);
        end
        [xcorrs, delays] = masked_norm_xcorr_no_fft_specific_delays(a, b, aBitmask,  bBitmask, aIsCircular, bIsCircular, [minDelay, maxDelay], noFFT);
        return;
    elseif (lenA < lenB) % a is not circular and is shorter
        a = zero_extend_linear_seq(a, 0, maxLen - lenA);
        aBitmask = zero_extend_linear_seq(aBitmask, 0, maxLen - lenA);
    end
    % by this point either an error should have been thrown
    %  or a and b should be the same length, n,
    %  and if either is circular, b is circular


   [xcorrs, delays, coverageLens] = unvalidated_masked_norm_xcorr(a, b, aBitmask, bBitmask, bIsCircular, minDelay, maxDelay, noFFT);

    if asRow
        xcorrs = xcorrs';
        coverageLens = coverageLens';
    end
end

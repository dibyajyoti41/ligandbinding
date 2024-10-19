function [sdps, delays] = lin_sdps_fft(linSeqA, linSeqB, minDelay, maxDelay)
    % lin_sdps_fft - Computes regular (linear) sliding dot products
    %  (a.k.a. raw cross-correlation coefficients) within a delay
    %  range using fourier transforms
    %
    % Inputs:
    %  linSeqA & linSeqB (column vectors representing the
    %    sequences to compute the cross-correlation over
    %    neither sequence is repeating -- treat them as if they
    %    are indefinitely preceded and followed by zeros)
    %  minDelay (non-positive value representing the minimum delay
    %    we are interested in for the range of delays)
    %  maxDelay (non-negative value representing the maximum delay
    %    we are interested in for the range of delays)
    % Outputs:
    %  sdps (column vector containing raw cross-correlations
    %    -- a.k.a. sliding dot products -- for the delay range)
    %  delays (row vector containing the delay which
    %    corresponds to each index in sdps)
    import SignalRegistration.Xcorr.zero_extend_linear_seq;

    validateattributes(linSeqA, {'numeric', 'logical'}, {'column', 'nonempty', 'real', 'finite', 'nonnan'});
    validateattributes(linSeqB, {'numeric', 'logical'}, {'column', 'nonempty', 'real', 'finite', 'nonnan'});

    lenA = length(linSeqA);
    lenB = length(linSeqB);
    n = max(lenA, lenB);
    if nargin < 3
        minDelay = 1 - n;
    else
        validateattributes(minDelay, {'numeric'}, {'scalar', 'integer', '<=', 0});
    end
    if nargin < 4
        maxDelay = n - 1;
    else
        validateattributes(maxDelay, {'numeric'}, {'scalar', 'integer', 'nonnegative'});
    end
    maxAbsDelay = max(abs([minDelay, maxDelay]));
    pow2 = 2^nextpow2(2*n - 1); % supposedly zero-padding sequence lengths to a power of 2 (at the end) can improve speed
    linSeqA = zero_extend_linear_seq(linSeqA, maxAbsDelay, pow2 - maxAbsDelay - lenA);
    linSeqB = zero_extend_linear_seq(linSeqB, maxAbsDelay, pow2 - maxAbsDelay - lenB);
    sdps = real(ifft(fft(linSeqA, pow2).*conj(fft(linSeqB, pow2))));
    lenDiff = maxAbsDelay - n;
    if lenDiff >= 0
        sdps = zero_extend_linear_seq([sdps(end-n+2:end); sdps(1:n)], lenDiff + 1, lenDiff + 1);
    else
        sdps = [sdps((end-maxAbsDelay+1):end); sdps(1:(maxAbsDelay+1))];
    end
    delays = minDelay:maxDelay;
    sdps = sdps(delays + ((end + 1)/2));
end

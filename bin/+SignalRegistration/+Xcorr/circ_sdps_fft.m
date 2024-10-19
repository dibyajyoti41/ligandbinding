function [sdps, delays] = circ_sdps_fft(refSeq, circSeqPeriod, minDelay, maxDelay)
    % circ_sdps_fft - Computes "circular" sliding dot products
    %  (a.k.a. raw cross-correlation coefficients) within a delay
    %  range using fourier transforms
    %  
    % Inputs:
    %  refSeq & circSeqPeriod (column vectors representing the
    %    sequences to compute the cross-correlation over --
    %    circSeqPeriod must represent a period from a repeating
    %    sequence since it will be treated as if it is indefinitely
    %    preceded and followed by it's own sequence)
    %  minDelay (non-positive value representing the minimum delay
    %    we are interested in for the range of delays)
    %  maxDelay (non-negative value representing the maximum delay
    %    we are interested in for the range of delays)
    %
    % Outputs:
    %  sdps (column vector containing raw cross-correlations
    %    -- a.k.a. sliding dot products -- for the delay range)
    %  delays (row vector containing the delay which
    %    corresponds to each index in sdps)
    import SignalRegistration.Xcorr.zero_extend_linear_seq;
    import SignalRegistration.Xcorr.extend_circular_seq;

    validateattributes(refSeq, {'numeric', 'logical'}, {'column', 'nonempty', 'real', 'finite', 'nonnan'});
    validateattributes(circSeqPeriod, {'numeric', 'logical'}, {'column', 'nonempty', 'real', 'finite', 'nonnan'});

    refLen = length(refSeq);
    circLen = length(circSeqPeriod);
    if (circLen < refLen)
        error('Cannot compute cross-correlations when the circular sequence is of shorter length than the reference sequence.');
    end
    if nargin < 3
        minDelay = 1 - circLen;
    else
        validateattributes(minDelay, {'numeric'}, {'scalar', 'integer', '<=', 0});
    end
    if nargin < 4
        maxDelay = circLen - 1;
    else
        validateattributes(maxDelay, {'numeric'}, {'scalar', 'integer', 'nonnegative'});
    end
    maxAbsDelay = max(abs([minDelay, maxDelay]));

    % should pad by maxAbsDelay because we want delays to take into account the periodic nature of circSeq
    % supposedly padding sequence lengths to a power of 2 can improve performance
    pow2 = 2^nextpow2(circLen + 2*maxAbsDelay);
    refSeq = zero_extend_linear_seq(refSeq, maxAbsDelay, pow2 - maxAbsDelay - refLen);
    circSeqPeriod = extend_circular_seq(circSeqPeriod, maxAbsDelay, pow2 - maxAbsDelay - circLen);

    sdps = real(ifft(fft(refSeq).*conj(fft(circSeqPeriod))));
    lenDiff = maxAbsDelay - circLen;
    % get rid of junk values from zero-padded computations and
    % correct delay ordering
    if lenDiff >= 0
        sdps = extend_circular_seq(sdps(1:circLen), maxAbsDelay, lenDiff + 1);
    else
        sdps = [sdps((end-maxAbsDelay+1):end); sdps(1:(maxAbsDelay+1))];
    end
    delays = minDelay:maxDelay;
    sdps = sdps(delays + ((end + 1)/2));
end
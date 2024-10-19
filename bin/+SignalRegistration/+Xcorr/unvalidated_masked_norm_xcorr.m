function [xcorrs, delays, coverageLens] = unvalidated_masked_norm_xcorr(a, b, aBitmask, bBitmask, bIsCircular, minDelay, maxDelay, noFFT)
    % UNVALIDATED_MASKED_NORM_XCORR - Computes pearson cross-correlations
    %   coefficients within a delay range
    %
    % SEE MASKED_NORM_XCORR
    %
    % Note: this works the same as MASKED_NORM_XCORR
    %  except it makes the assumption that all the inputs are
    %  provided and that they are all "valid":
    %  a & b must be nonempty column vectors of real,
    %   finite, nonnan numbers
    %  aBitmask & bBitmask must be logical vectors of the same
    %   size as a & b
    %  bIsCircular must be true or false and if it is false,
    %   neither a nor b ought to be circular and if it is true,
    %   b must not be shorter than a
    %  minDelay must be a non-positive integer
    %  maxDelay must be a non-negative integer
    %  noFFT must be true or false
    %  In cases where performance is critical this function can be
    %  used after making sure that all the inputs provided will 
    %  always be valid.
    import SignalRegistration.Xcorr.masked_norm_xcorr_fft;
    import SignalRegistration.Xcorr.masked_norm_xcorr_no_fft_specific_delays;
% 
%     tmpStruct = struct();
%     tmpStruct.a = a;
%     tmpStruct.b = b;
%     tmpStruct.aBitmask = aBitmask;
%     tmpStruct.bBitmask = bBitmask;
%     tmpStruct.bIsCircular = bIsCircular;
%     tmpStruct.minDelay = minDelay;
%     tmpStruct.maxDelay = maxDelay;
%     tmpStruct.noFFT = noFFT;
%     

    DEFAULT_FFT_THRESHOLD = 0; %Todo: optimize (currently always uses it unless noFFT is true)
    if noFFT
        FFT_THRESHOLD = Inf;
    else
        FFT_THRESHOLD = DEFAULT_FFT_THRESHOLD;
    end

    if length(b) > FFT_THRESHOLD
        [xcorrs, delays, coverageLens] = masked_norm_xcorr_fft(a, b, aBitmask, bBitmask, bIsCircular, minDelay, maxDelay);
    else
        [xcorrs, delays, coverageLens] = masked_norm_xcorr_no_fft_specific_delays(a, b, aBitmask, bBitmask, bIsCircular, minDelay:maxDelay);
    end
% 
%     tmpStruct.xcorrs = xcorrs;
%     
%     tmpStruct.delays = delays;
%     tmpStruct.coverageLens = coverageLens;
% 
%     assignin('base', 'tmpStructV2', tmpStruct);
end

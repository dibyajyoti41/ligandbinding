function [maxPossibleCoverageLen] = get_max_possible_coverage_len(bitmaskA, bitmaskB, circularA, circularB, delay)
    % get_max_possible_coverage_len
    %  Computes the maximum coverage length one can achieve by
    %  aligning sequences with the provided bitmasks
    %  up to a given delay
    %
    % Inputs:
    %  bitmaskA & bitmaskB - bitmasks for sequences
    %  circularA & circularB - whether the sequences
    %   given should be treated like they are circular
    %   (e.g. allow circularly shifting them, respecting their
    %     periodic nature) as opposed to linear sequences
    %    (which cannot be circularly shifted and don't repeat)
    %  delay - the magnitude of delay (shifting/offsetting)
    %    permitted
    %
    % Outputs:
    %  maxPossibleCoverageLen - maximum coverage length possible
    %    for the bitmasks in question
    import SignalRegistration.Xcorr.circ_sdps_fft;
    import SignalRegistration.Xcorr.lin_sdps_fft;

    maxLen = max(length(bitmaskA), length(bitmaskB));
    if nargin < 5
        delay = maxLen - 1;
    else
        delay = abs(delay);
    end
    delay = [-1*delay, delay];
    minDelay = delay(1);
    maxDelay = delay(2);

    if (circularA && not(circularB))
        tmpBitmask = bitmaskA;
        bitmaskA = bitmaskB;
        bitmaskB = tmpBitmask;
        % circularA = false;
        circularB = true;
    end
    bitmaskA = bitmaskA(:);
    bitmaskB = bitmaskB(:);

    if circularB
        lenA = length(bitmaskA);
        lenB = length(bitmaskB);
        minLen = min(lenA, lenB);
        lenAStarts = 1;
        if (lenA > minLen)
            lenAStarts = 1:(1 + lenA - minLen);
        end
        maxPossibleCoverageLen = 0;
        for lenAStart = lenAStarts
            bitmaskACropped = bitmaskA(lenAStart:(lenAStart + minLen - 1));
            maxPossibleCoverageLen = max(maxPossibleCoverageLen, round(max(circ_sdps_fft(bitmaskACropped, bitmaskB, minDelay, maxDelay))));
        end
    else
        maxPossibleCoverageLen = round(max(lin_sdps_fft(bitmaskA, bitmaskB, minDelay, maxDelay)));
    end
end
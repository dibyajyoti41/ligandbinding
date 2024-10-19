function [scores] = compute_scores(xcorrs, coverageLens, maxPossibleCoverageLen) %#ok<INUSD>
    % compute_scores
    %  Computes the scores array based off of the xcorrs,
    %   coverageLens, and maxCoverageLen
    %
    % Inputs:
    %  xcorrs - raw Pearson cross-correlation coefficient
    %    values as computed in a three-dimensional array.
    %    The first dimensions specifies whether a sequence was
    %     flipped (index 2 if true, 1 if false).
    %    The second dimension specifies how much circular
    %      shifting there was (index 1 = no shift, 2 = shift of 1,
    %    etc.)
    %    The third dimension specifies the linear offset,
    %      (index k => firstOffset + k linear shifting)
    %  coverageLens - Effective sample sizes corresponding to each
    %    xcorrs value
    %  maxPossibleCoverageLen - maximum coverage length possible
    %    for the bitmasks in question
    %
    % Outputs:
    %  scores - Weighted scores for xcorrs/coverageLens
    %    (a function with a direct relationship to
    %      Pearson cross-correlation coefficient
    %      and the coverage length in proportion to the
    %      maximum coverage length possible for the bitmasks
    %      in question) 
    %

    scores = xcorrs(:).*sqrt(coverageLens(:));
    scores = reshape(scores, size(xcorrs));
end
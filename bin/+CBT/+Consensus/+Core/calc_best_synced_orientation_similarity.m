function [bsosStruct] = calc_best_synced_orientation_similarity(circSeqA, circSeqB, indexWeightsA, indexWeightsB)
    import SignalRegistration.XcorrAlign.get_no_crop_lin_circ_xcorrs;
    import SignalRegistration.XcorrAlign.compute_scores;
    import SignalRegistration.XcorrAlign.get_max_possible_coverage_len;

    if (length(circSeqA) ~= length(circSeqB))
        error('Sequences should be the same length!');
    end

    bitmaskA = logical(indexWeightsA);
    bitmaskB = logical(indexWeightsB);
    [xcorrs, coverageLens] = get_no_crop_lin_circ_xcorrs(circSeqA, circSeqB, bitmaskA, bitmaskB);
    scores = compute_scores(xcorrs, coverageLens, get_max_possible_coverage_len(bitmaskA, bitmaskB, true, true));


    % we define the best synced orientation (BSO) as the reorientation
    % that produces the highest score based on the pearson
    % correlation coefficient and the coverage (sample size) of the
    % sequences
    %
    % in this case both sequences are assumed to be cyclical and
    %  containing comparable numbers of basepairs per
    %  datapoint/pixel
    % all pairs of cyclical permutations for the flipped and
    %  nonflipped versions of the sequences are explored
    % coverage can vary if the index weights vary
    [bestScore, bestScoreLinearIndex] = max(scores(:));
    [idxFlipping, idxCircshifting, ~] = ind2sub(size(scores), bestScoreLinearIndex);
    xcorrAtBest = xcorrs(idxFlipping, idxCircshifting, 1);

    flipTFAtBest = logical(idxFlipping - 1);
    circShiftAtBest = 0 - (idxCircshifting - 1);
    % if B is circshifted with the value of circShiftAtBest
    %  and then flipped (if flipTFAtBest is true)
    % it should look synced to A


    bsosStruct = struct();

    % Store similarity scoring data
    bsosStruct.bestScore = bestScore; % orientation was synced to maximize this similarity score
    bsosStruct.xcorrAtBest = xcorrAtBest; % pearson correlation coefficient associated with the synced orientation

    % Store relative orientation data
    bsosStruct.flipTFAtBest = flipTFAtBest;
    bsosStruct.circShiftAtBest = circShiftAtBest;
end
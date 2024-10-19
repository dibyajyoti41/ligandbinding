function [signalStartIdx, signalEndIdx] = find_signal_region_naive(curve)
    % FIND_SIGNAL_REGION_NAIVE - finds the "signal region" of a molecule intensity 
    %	trace (i.e., the part of the curve that corresponds to a molecule and
    %	not to background)
    %
    % Inputs:
    %   curve
    %     the curve
    %
    % Outputs:
    %   signalStartIdx
    %     the index of the approximated left boundary for the molecule
    %   signalEndIdx
    %     the index of the approximated right boundary for the molecule
    %
    % Authors:
    %   Charleston Noble

    curve = abs(curve);

    mCurve = nanmean(curve);

    signalStartIdx = max([1, find(curve > mCurve, 1, 'first')]);
    signalEndIdx = min([length(curve) find(curve > mCurve, 1, 'last')]);
end
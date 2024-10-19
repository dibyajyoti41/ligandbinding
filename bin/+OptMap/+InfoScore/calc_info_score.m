function [infoScore] = calc_info_score(energyLandscape, energyBarrierThreshold, estNoiseStd, smoothingWindowSz)
    % CALC_INFO_SCORE - calculates the information score for an
    %	 energyLandscape given an energyBarrierThreshold
    %  See the WPAlign paper for details
    %
    % Inputs:
    %   energyLandscape
    %     the energy landscape curve
    %   energyBarrierThreshold
    %     the minimum difference in values between adjacent robust
    %     extrema in the landscape to be used for calculating the
    %     information score
    %   estNoiseStd
    %     the standard deviation of the estimated noise
    %   smoothingWindowSz
    %
    % Outputs:
    %   infoScore
    %     the information score
    %
    % Dependencies:
    %   OptMap.SignalProcessing.detect_robust_local_extrema
    %
    % Authors:
    %  Saair Quaderi
    import OptMap.SignalProcessing.detect_robust_local_extrema;

    if nargin < 4
        smoothingWindowSz = 1;
    else
        validateattributes(smoothingWindowSz, {'numeric'}, {'scalar', 'positive', 'integer'});
    end

    if smoothingWindowSz > 1
        % Smooth the curve
        energyLandscape = smooth(energyLandscape, smoothingWindowSz);
    end

    % SI suggests a different approach (a modified version of Azbel's
    %   method) to getting the robust local extrema, but I think this
    %   approach should work at least as effectively
    [~, extremaIntensityVals] = detect_robust_local_extrema(energyLandscape, energyBarrierThreshold);

    % Get the diffs of the extrema
    robustIntensityDiffs = abs(diff(extremaIntensityVals));

    % Calculate information
    x = log(robustIntensityDiffs);
    xSquaredDistFromMu = x.^2; % (x - mu).^2 where mu=0

    chiRegularizationParam = 1;

    logOfSigmaSquaredPlusChi = log(estNoiseStd.^2 + chiRegularizationParam);

    probVals = (1 ./ (sqrt(2*pi.*logOfSigmaSquaredPlusChi))) .* exp( -xSquaredDistFromMu./(2 .* logOfSigmaSquaredPlusChi));
    % functionally equivalent to y = normpdf(x, mu, sigma);
    %  with mu = 0 and sigma=sqrt(sigmaSquared)

    probVals = probVals(probVals > 0);
    infoScore = sum(-log(probVals));
end
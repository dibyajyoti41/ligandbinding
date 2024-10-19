function [kymoInfoScore] = calc_kymo_info_score(alignedKymoFgImg, numSigma, smoothingWindowSz)
    % CALC_KYMO_INFO_SCORE - calculates the information
    %	 score for an intensity curve calculated from an aligned
    %    kymograph
    %  See the WPAlign paper for details
    %
    % Inputs:
    %   alignedKymoFgImg
    %     the foreground portion of the aligned kymograph
    %   numSigma
    %     (optional, defaults to 1)
    %     multiple of sigma that is to be used as an energy barrier
    %     threshold when detecting extrema for info score
    %   smoothingWindowSz
    %
    % Outputs:
    %   kymoInfoScore
    %     the information score for the kymograph
    %
    % Authors:
    %   Saair Quaderi
    import OptMap.InfoScore.calc_info_score;

    validateattributes(alignedKymoFgImg, {'numeric'}, {'nonempty', '2d', 'real', 'finite', 'nonnegative'});

    if nargin < 2
        numSigma = 1;
    else
        validateattributes(numSigma, {'numeric'}, {'scalar', 'real', 'finite', 'nonnegative'});
    end

    if nargin < 3
        smoothingWindowSz = 1;
    else
        validateattributes(smoothingWindowSz, {'numeric'}, {'scalar', 'positive', 'integer'});
    end

    % Smooth the curve
    energyLandscape =  mean(alignedKymoFgImg);
    if (smoothingWindowSz > 1)
         % Note: Paper/SI doesn't specify smoothing
        energyLandscape = smooth(energyLandscape, smoothingWindowSz);
    end

    % Calculate the threshold
    estSignalImgNoise = alignedKymoFgImg - repmat(energyLandscape, [size(alignedKymoFgImg, 1), 1]);

    %  SI suggests a value of <(estImgNoise - <estImgNoise>)^2>
    %    for sigma squared and calls it the background variance
    % not quite sure whether there is a reason why normalizing by
    %    N (the sample size) would actually be preferable
    %    over normalizing by N-1 for the variance, but that's why
    %    we're using std(estSignalImgNoise(:), 1), the square root
    %     of the second moment of the sample about its mean
    sigma_estSignalImgNoiseStd = std(estSignalImgNoise(:), 1);

    % SI says "typically chosen equal to the average background noise, i.e., sigma"
    %  which seems to be a contradiction since sigma is not an estimate
    %  of the average noise but rather one of the std dev of the noise
    %  some multiple of std dev seems more sensible so choosing the value below:
    energyBarrierThreshold = numSigma*sigma_estSignalImgNoiseStd;

    % Note: info score calculation uses a somewhat different method
    %   for detecting robust local extrema than the SI suggests
    kymoInfoScore = calc_info_score(energyLandscape, energyBarrierThreshold, sigma_estSignalImgNoiseStd, 1);
end
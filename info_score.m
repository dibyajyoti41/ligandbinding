 function [infoScore] = info_score(noiseThreshold,xMinima,xMaxima)

    xMinima = xMinima(2:end);
     % Get the diffs of the extrema
    robustIntensityDiffs = abs(diff(sort([xMinima; xMaxima ])));

    % Calculate information
    x = log(robustIntensityDiffs);
    xSquaredDistFromMu = x.^2; % (x - mu).^2 where mu=0

    chiRegularizationParam = 1;


    logOfSigmaSquaredPlusChi = log(noiseThreshold.^2 + chiRegularizationParam);

    probVals = (1 ./ (sqrt(2*pi.*logOfSigmaSquaredPlusChi))) .* exp( -xSquaredDistFromMu./(2 .* logOfSigmaSquaredPlusChi));
  
    % functionally equivalent to y = normpdf(x, mu, sigma);
    %  with mu = 0 and sigma=sqrt(sigmaSquared)

%     probVals = probVals(probVals > 0);
    infoScore = sum(-log(probVals));
% 

end


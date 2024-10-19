function [fittingModel, fittingStartPoints, fittingLowerBounds, fittingUpperBounds] = get_double_tanh_fitting_model_info(minIdx, maxIdx, minVal, maxVal, minTanhStretchFactor, maxTanhStretchFactor, tanhStretchFactorInitGuess, startEdgeIdxInitGuess, endEdgeIdxInitGuess, heightDiffInitGuess)
    % GET_DOUBLE_TANH_FITTING_MODEL_INFO - supplies appropriate
    %   fitting model information, including constraints and
    %   initial approximations using the following function
    %   f(x) = A + F *(tanh( (x - B) * D ) - tanh( (x - C) * E ))
    %
    % Note:
    %  Be careful to deal with unintended fits where B >= C
    %    since that isn't specifically constrained against
    %
    %
    % Inputs:
    %   minIdx
    %   maxIdx
    %   minVal
    %   maxVal
    %   minTanhStretchFactor
    %   maxTanhStretchFactor
    %   tanhStretchFactorInitGuess
    %   startEdgeIdxInitGuess
    %   endEdgeIdxInitGuess
    %   heightDiffInitGuess
    %
    % Outputs:
    %   fittingModel
    %   fittingStartPoints
    %   fittingLowerBounds
    %   fittingUpperBounds
    %
    % Authors:
    %   Saair Quaderi
    %    (complete refactoring & adjusting some constraints)
    %   Charleston Noble
    %    (earlier version)

    if nargin < 8
        startEdgeIdxInitGuess = minIdx;
    end
    if nargin < 9
        endEdgeIdxInitGuess = maxIdx;
    end
    if nargin < 10
        heightDiffInitGuess = maxVal - minVal;
    end

    fittingModel = fittype('A + F * (tanh( (x - B) * D ) - tanh( (x - C) * E ))');


    % --- Define lower and upper bounds to constrain the fit ---

    % The lower limit for values that can be produced by the 
    %  fitting function
    lowerBoundA = minVal; % absolute lower bound
    upperBoundA = maxVal; % absolute upper bound, but actual value should be lower in practice

    % Horizontal central index position of the left edge, tanh
    lowerBoundB = minIdx; % absolute lower bound, assuming left edge is in view
    upperBoundB = maxIdx; % absolute upper bound, assuming left edge is in view

    % Horizontal central index position of the right edge, -tanh
    lowerBoundC = minIdx; % absolute lower bound, assuming right edge is in view
    upperBoundC = maxIdx; % absolute upper bound, assuming right edge is in view

    % Factor for horizontal stretching of the left edge, tanh
    %  (smaller value means the curve rises more steeply)
    lowerBoundD = minTanhStretchFactor;
    upperBoundD = maxTanhStretchFactor;

    % Factor for horizontal stretching of the right edge, -tanh
    %  (smaller value means the curve falls more steeply)
    lowerBoundE = minTanhStretchFactor;
    upperBoundE = maxTanhStretchFactor;

    % Half the distance between the lower and upper limit 
    %  for values that can be produced by the fitting
    %  function
    lowerBoundF = 0.0; % absolute lower bound, but actual value should be higher in practice
    upperBoundF = (maxVal - minVal)/2.0; % absolute upper bound

    % -- Set initial guesses for values --
    initGuessA = lowerBoundA;
    initGuessB = startEdgeIdxInitGuess;
    initGuessC = endEdgeIdxInitGuess;
    initGuessD = tanhStretchFactorInitGuess;
    initGuessE = tanhStretchFactorInitGuess;
    initGuessF = heightDiffInitGuess/2.0;

    fittingStartPoints = [initGuessA, initGuessB, initGuessC, initGuessD, initGuessE, initGuessF];
    fittingLowerBounds = [lowerBoundA, lowerBoundB, lowerBoundC, lowerBoundD, lowerBoundE, lowerBoundF];
    fittingUpperBounds = [upperBoundA, upperBoundB, upperBoundC, upperBoundD, upperBoundE, upperBoundF];
end
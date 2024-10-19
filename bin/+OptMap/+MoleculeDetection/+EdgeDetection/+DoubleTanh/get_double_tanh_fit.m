function [fittingResult, gofStruct, outputStruct] = get_double_tanh_fit(valVect, moleculeStartEdgeIdxFirstApprox, moleculeEndEdgeIdxFirstApprox, tanhSettingsStruct)
    % GET_DOUBLE_TANH_FIT - Fits the curve against the following
    %    function, using the specified locations for the start
    %    and end edges (associated with B and C) as the initial
    %    guesses
    %  f(x) = A + F *(tanh( (x - B) * D ) - tanh( (x - C) * E ))
    %
    % Inputs:
    %   valVect
    %   moleculeStartEdgeIdxFirstApprox
    %   moleculeEndEdgeIdxFirstApprox
    %   tanhSettingsStruct
    %
    % Outputs:
    %   fittingResult
    %   gofStruct
    %   outputStruct
    %
    % Authors:
    %   Saair Quaderi
    %   Charleston Noble
    %     (earlier version)

    import OptMap.MoleculeDetection.EdgeDetection.DoubleTanh.get_double_tanh_fitting_model_info;

    vectLen = length(valVect);
    try
    minNrmMoleculeValApprox = min(valVect(moleculeStartEdgeIdxFirstApprox:moleculeEndEdgeIdxFirstApprox));
    maxNrmMoleculeValApprox = max(valVect(moleculeStartEdgeIdxFirstApprox:moleculeEndEdgeIdxFirstApprox));
    catch e
        error(e)
    end
    [fittingModel, fittingStartPoints, fittingLowerBounds, fittingUpperBounds] = get_double_tanh_fitting_model_info(...
        1, vectLen,...
        minNrmMoleculeValApprox, maxNrmMoleculeValApprox,...
        tanhSettingsStruct.minTanhStretchFactor, tanhSettingsStruct.maxTanhStretchFactor,...
        tanhSettingsStruct.tanhStretchFactorInitGuess,...
        moleculeStartEdgeIdxFirstApprox, moleculeEndEdgeIdxFirstApprox);

    % Fit to the model
    [fittingResult, gofStruct, outputStruct] = fit( ...
        (1:vectLen)', valVect(:), ...
        fittingModel, ...
        'StartPoint', fittingStartPoints, ...
        'Lower', fittingLowerBounds, ...
        'Upper', fittingUpperBounds);
end
function [unroundedMoleculeStartEdgeIdxs, unroundedMoleculeEndEdgeIdxs] = adjust_nrm_kymo_edge_detection_unrounded(kymoFramewiseNrm, moleculeStartEdgeIdxsFirstApprox, moleculeEndEdgeIdxsFirstApprox, tanhSettings)
    % ADJUST_NRM_KYMO_EDGE_DETECTION_UNROUNDED - attempts
    %  to improve edge detection from given initial approximations
    %  by fitting the intensity profiles associated with each frame in
    %  the kymograph to the following function and extracting the
    %  B and C positions as the start and end edge coordinates:
    %  f(x) = A + F *(tanh( (x - B) * D ) - tanh( (x - C) * E ))
    %
    % Inputs:
    %   kymoFramewiseNrm
    %   moleculeStartEdgeIdxsFirstApprox
    %   moleculeEndEdgeIdxsFirstApprox
    %   tanhSettings
    %
    % Outputs:
    %   unroundedMoleculeStartEdgeIdxs
    %   unroundedMoleculeEndEdgeIdxs
    %
    % Authors:
    %   Saair Quaderi
    %

    import OptMap.MoleculeDetection.EdgeDetection.DoubleTanh.get_double_tanh_fit;

    numFrames = size(kymoFramewiseNrm, 1);
    frameNums = (1:numFrames)';
    fittingResults = arrayfun(@(frameNum) ...
        get_double_tanh_fit(...
            kymoFramewiseNrm(frameNum, :),...
            moleculeStartEdgeIdxsFirstApprox(frameNum),...
            moleculeEndEdgeIdxsFirstApprox(frameNum),...
            tanhSettings),...
        frameNums, ...
        'UniformOutput', false);
    unroundedMoleculeStartEdgeIdxs = cellfun(@(fittingResult) fittingResult.B, fittingResults);
    unroundedMoleculeEndEdgeIdxs = cellfun(@(fittingResult) fittingResult.C, fittingResults);

    badFits = unroundedMoleculeStartEdgeIdxs >= unroundedMoleculeEndEdgeIdxs;
    if any(badFits)
        warning('Some invalid double-tanh fits were detected and reverted back to first approximations');
        unroundedMoleculeStartEdgeIdxs(badFits) = moleculeStartEdgeIdxsFirstApprox(badFits);
        unroundedMoleculeEndEdgeIdxs(badFits) = moleculeEndEdgeIdxsFirstApprox(badFits);
    end
end
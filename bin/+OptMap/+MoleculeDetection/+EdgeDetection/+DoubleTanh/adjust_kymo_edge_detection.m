function [moleculeStartEdgeIdxs, moleculeEndEdgeIdxs, mainKymoMoleculeMask] = adjust_kymo_edge_detection(kymo, moleculeStartEdgeIdxsFirstApprox, moleculeEndEdgeIdxsFirstApprox, tanhSettings)
     % ADJUST_KYMO_EDGE_DETECTION_UNROUNDED - attempts
    %  to improve edge detection from given initial approximations
    %  by fitting the intensity profiles associated with each frame in
    %  the kymograph to the following function and extracting the
    %  B and C positions as the start and end edge coordinates:
    %  f(x) = A + F *(tanh( (x - B) * D ) - tanh( (x - C) * E ))
    %
    % Inputs:
    %   kymo
    %   moleculeStartEdgeIdxsFirstApprox
    %   moleculeEndEdgeIdxsFirstApprox
    %   tanhSettings
    %
    % Outputs:
    %   moleculeStartEdgeIdxs
    %   moleculeEndEdgeIdxs
    %
    % Authors:
    %   Saair Quaderi
    %

    kymoFrameMinVals = min(kymo, [], 2);
    kymoFrameMaxVals = max(kymo, [], 2);
    kymoFramewiseNrm = bsxfun(@rdivide, bsxfun(@minus, kymo, kymoFrameMinVals), (kymoFrameMaxVals - kymoFrameMinVals));
    kymoFramewiseNrm(isnan(kymoFramewiseNrm)) = 0;
    
    import OptMap.MoleculeDetection.EdgeDetection.DoubleTanh.adjust_nrm_kymo_edge_detection_unrounded;
    [unroundedMoleculeStartEdgeIdxs, unroundedMoleculeEndEdgeIdxs] = adjust_nrm_kymo_edge_detection_unrounded(...
        kymoFramewiseNrm, ...
        moleculeStartEdgeIdxsFirstApprox, moleculeEndEdgeIdxsFirstApprox, ...
        tanhSettings);
    moleculeStartEdgeIdxs = round(unroundedMoleculeStartEdgeIdxs);
    moleculeEndEdgeIdxs = round(unroundedMoleculeEndEdgeIdxs);
    mainKymoMoleculeMask = false(size(kymo));
    numFrames = size(kymo, 1);
    for frameNum = 1:numFrames
        mainKymoMoleculeMask(frameNum, moleculeStartEdgeIdxs(frameNum):moleculeEndEdgeIdxs(frameNum)) = true;
    end
end
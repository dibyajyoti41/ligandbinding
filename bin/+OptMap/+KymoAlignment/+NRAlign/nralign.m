function [alignedKymo, stretchFactorsMat, shiftAlignedKymo, alignedMask, shiftAlignedKymoMask, shiftingVect, alignmentSuccessTF] = nralign(unalignedKymo, skipPrealignTF, kymoMask)
    % NRALIGN - non-recursively aligns a kymograph, using a modified
    %  wpalign algorithm.
    %
    % Inputs:
    %	unalignedKymo
    %	  the unaligned kymograph
    %   skipPrealignTF
    %     true if pre-alignment step should be skipped (e.g. if it has
    %     already been pre-aligned)
    %   kymoMask
    %     true wherever the foreground element of interest is believed to
    %     be located (helps ignore values in background that aren't
    %     relevant)
    %
    % Outputs:
    %	alignedKymo
    %	  the aligned kymograph is output in the last step
    %   stretchFactorsMat
    %	  array of the stretch factor of each pixel in the kymograph
    %   shiftAlignedKymo
    %     the shift aligned kymograph on which stretch factors were
    %     computed to produce the aligned kymograph
    %   alignedMask
    %     the aligned mask
    %   shiftAlignedKymoMask
    %     the shift aligned mask
    %   shiftingVect
    %     the amount each row was shifted prior to feature-based alignment
    %   alignmentSuccessTF
    %     whether the alignment attempt was succesful, false if features
    %     could not be found
    %
    % Authors:
    %	Henrik Nordanger
    %   Saair Quaderi

    if nargin < 2
        skipPrealignTF = false;
    end
    if nargin < 3
        kymoMask = true(size(unalignedKymo));
    end
    
    % %   forceEdgesTF (optional; defaults to false)
    % %     true if forcing alignment near edges
    % if nargin < 4
    %     forceEdgesTF = false;
    % end
    % edgeFeatureAlignmentBuffer = 10;
    
    shiftAlignedKymo = unalignedKymo;
    shiftAlignedKymoMask = kymoMask;
    if not(skipPrealignTF)
        import OptMap.KymoAlignment.NRAlign.pre_nralign_shift_align;
        [shiftAlignedKymo, shiftingVect] = pre_nralign_shift_align(shiftAlignedKymo);
        import OptMap.KymoAlignment.NRAlign.shift_rows;
        shiftAlignedKymoMask = shift_rows(shiftAlignedKymoMask, shiftingVect);
    else
        shiftingVect = zeros(size(shiftAlignedKymo, 1));
    end
    
    
    % The number of pixels around each feature on each side (not including
    % the pixel the feature is on) on which there may not be another
    % feature
    W_trim = 5;
    typicalFeatureWidth = (2 * W_trim) + 1;

    % The maximum number of pixels a feature can move from frame to next
    % frame
    maxFeatureMovementPx = 3;


    % The number of pixels to ignore at the edges of the kymograph
    %  when looking for features
    featurelessSideWidth = 0;
    % The edges of the region in which features are to be found are
    %  specified, unless already done in the call to the function.
    leftColIdx = 1 + featurelessSideWidth;
    rightColIdx = size(unalignedKymo, 2) - featurelessSideWidth;


    % The (maximum) number of features to look for "k"
    minFgMaskWidth = min(min(sum(shiftAlignedKymoMask, 2)), (size(shiftAlignedKymoMask, 2) - 2*featurelessSideWidth));
    maxNumFeaturesSoughtK = ceil(minFgMaskWidth / ((2 * W_trim) + 1));

    
    % The shift-aligned kymograph is smoothed out in order to reduce noise
    squareSmoothingWindowLen_pixels = 10;
    blurSigmaWidth_pixels = 2;
    
    import OptMap.KymoAlignment.apply_gaussian_blur;
    smoothImg = apply_gaussian_blur(shiftAlignedKymo, squareSmoothingWindowLen_pixels, blurSigmaWidth_pixels);


    %The laplacian of gaussian is obtained and normalized
    import OptMap.KymoAlignment.apply_laplacian_of_gaussian_filter;
    smoothImg = apply_laplacian_of_gaussian_filter(smoothImg, [2, 6], 2);
    smoothImg = smoothImg ./ max(abs(smoothImg(:)));

    fgMaskExpansionSideWidth = 10;
    if not(all(shiftAlignedKymoMask(:)))
        farBgMask = ~imdilate(shiftAlignedKymoMask, true(1, 2*fgMaskExpansionSideWidth + 1));
        smoothImg(farBgMask) = NaN;
    end
    
    alignmentSuccessTF = all(sum(~isnan(smoothImg) & (smoothImg ~= 0), 2));
    stretchFactorsMat = ones(size(shiftAlignedKymo));
    alignedKymo = shiftAlignedKymo;
    alignedMask = shiftAlignedKymoMask;
    
    
    if alignmentSuccessTF
        %The k shortest paths through the (prealigned) kymograph are found
        import OptMap.KymoAlignment.find_k_features;
        pathsColIdxOffsets = find_k_features(smoothImg(:, leftColIdx:rightColIdx), maxFeatureMovementPx, typicalFeatureWidth, maxNumFeaturesSoughtK); 
        
        alignmentSuccessTF = not(isempty(pathsColIdxOffsets));
        if alignmentSuccessTF
            pathsColIdxOffsets = horzcat(pathsColIdxOffsets{:, 1});
            pathsColIdxs = leftColIdx + pathsColIdxOffsets - 1;

            
            % % Force edge alignment -- doesn't currently seem to
            % %  behave quite as expected (todo: debug)
            % if forceEdgesTF
            %     numFrames = size(pathsColIdxs, 1);
            %     forcedLeftFeature = NaN(numFrames, 1);
            %     forcedRightFeature = NaN(numFrames, 1);
            %     for frameNum = 1:numFrames
            %         tmp = find(alignedMask(frameNum, :), 1, 'first');
            %         if isempty(tmp)
            %             tmp = NaN;
            %         else
            %             tmp = max(1, tmp - edgeFeatureAlignmentBuffer);
            %         end
            %         forcedLeftFeature(frameNum) = tmp;
            % 
            %         tmp = find(alignedMask(frameNum, :), 1, 'last');
            %         if isempty(tmp)
            %             tmp = NaN;
            %         else
            %             tmp = min(size(smoothImg, 2), tmp + edgeFeatureAlignmentBuffer);;
            %         end
            %         forcedRightFeature(frameNum) = tmp;
            %     end
            %     import OptMap.KymoAlignment.nearest_nonnan;
            %     forcedLeftFeature = nearest_nonnan(forcedLeftFeature, 1);
            %     forcedRightFeature = nearest_nonnan(forcedRightFeature, size(smoothImg, 2));
            %     pathsColIdxs = [forcedLeftFeature, pathsColIdxs, forcedRightFeature];
            % end
            
            imgSz = size(alignedKymo);
            import OptMap.KymoAlignment.compute_stretch_factors;
            stretchFactorsMat = compute_stretch_factors(pathsColIdxs, imgSz);

            import OptMap.KymoAlignment.apply_horizontal_stretching;
            alignedKymo = apply_horizontal_stretching(alignedKymo, stretchFactorsMat);
            alignedMask = round(apply_horizontal_stretching(double(alignedMask), stretchFactorsMat)) > 0;
        end
    end
    
    % % For debugging reasons, the paths can be drawn in the (unaligned)
    % % kymograph.
    % 
    % 
    % 
    % szImg = size(unalignedKymo);
    % pathLabelsMat = zeroes(szImg);
    % numPaths = size(pathsColIdxs, 2);
    % pathRowIdxs = 1:szImg(1);
    % for pathNum = 1:numPaths
    %     pathColIdxs = pathsColIdxs(:, pathNum);
    %     pathLabelsMat(sub2ind(pathLinIdxs, pathRowIdxs, pathColIdxs)) = pathNum;
    % end
    % 
    % hFig = figure();
    % hPanel = uipanel('Parent', hFig);
    % hAxis = axes('Parent', hPanel);
    % import OptMap.KymoAlignment.UI.plot_features_overlay;
    % plot_features_overlay(hAxis, unalignedKymo, pathLabelsMat);
end
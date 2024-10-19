function [spalignedKymo,spalignedBitmask,cutKymo,cutMaskF,f] = spalign(kymo,bitmask, minOverlap,maxShift, skipAlign, detPeaks)
    % spalign  shifted peaks alignment
    %
    %       Args:
    %           kymo - kymograph
    %           bitmask - where the signal is in the kymograph
    %           minOverlap - minimum overlap between two rows
    %           maxShift - max allowed shift
    %           skipAlign - whether to skip shift align (should skip for
    %           short barcodes, since then the local alignment not
    %           trustworthy)
    %
    %       Returns:
    %           spalignedKymo - kymo aligned using spalign
    %           spalignedBitmask - bitmask aligned using spalign
    %
    %
    
    % aligns features in a kymograph
    if nargin < 3
        minOverlap = 300;
        maxShift = 20;
    end

    %% STEP 1: shift align with minimum overlap (shouldn't be too small to
    % avoid erranoreously matching small regions)
    import OptMap.KymoAlignment.SPAlign.apply_stretching;
    import OptMap.KymoAlignment.SPAlign.compute_rescale_factors;
    if ~skipAlign
        [cutKymo,cutMaskF] = shift_align(kymo,bitmask,minOverlap,maxShift);
    else
        cutKymo = kymo;
        cutMaskF = bitmask;
    end
    
        alignedkymo = double(cutKymo);
        alignedkymo(~cutMaskF) = nan;
    
%                figure;
%  imshowpair(imresize(bitmask,[200 500]),imresize(kymo,[200 500]), 'ColorChannels','red-cyan'  )

    
    %% STEP 2: find peaks
%     https://github.com/BorgwardtLab/Topf
%     skipPeakDet = 0;
    % TDA to make peaks more prominent
    if detPeaks
        maxNumFeaturesSoughtK = 20;
        minPeakDistance = 11;
        gfiltSigma = 5;
        [vec,Ws,Ypk,Xpk,Wpk] = find_mat_peaks(alignedkymo, gfiltSigma,maxNumFeaturesSoughtK, minPeakDistance);
        try
            %% STEP 3: track peaks

            distPar = 5;
            [sourcesLong,pathsColIdxs,longFeatures,sources] = track_peaks(alignedkymo,vec,Ws,distPar);
        catch
            pathsColIdxs = [];
        end
    else
        pathsColIdxs = [];
    end
    
       % Two extra features
    leftF =  arrayfun(@(frameNum) find(cutMaskF(frameNum, :), 1, 'first'), 1:size(cutMaskF,1));
    rightF = arrayfun(@(frameNum) find(cutMaskF(frameNum, :), 1, 'last'), 1:size(cutMaskF,1));

    %The paths are sorted.
    pathsColIdxs = sortrows(pathsColIdxs');
    pathsColIdxs = [leftF'  pathsColIdxs'  rightF'];
%         pathsColIdxs =pathsColIdxs';

    imgSz = size(alignedkymo);
    % import DBM4.compute_rescale_factors;
    stretchFactorsMat = compute_rescale_factors(pathsColIdxs, imgSz);

    if size(pathsColIdxs,2)<=1
        stretchFactorsMat = ones(size(alignedkymo));
    end

% catch
%     stretchFactorsMat = ones(size(alignedkymo));
%             
% end


      
%        
%     figure
%     hold on
%         for s = 1:length(sources)
%             plot( sources{s}(:, 2),sources{s}(:, 1),'|','linewidth', 2)
%         end
%     sourcesLong
% img = alignedKymoF(:,st:stop);

    % alternativ: old
% alignedMask
%     [pathsColIdxs ] = align_old(alignedkymo,cutMaskF);
    

% first andlast features should not be stretched
% meanFirst = round(nanmean(pathsColIdxs(:,1)));
% for s = 1:size(pathsColIdxs,1)
%     
% end

 
    

% 
% stretchFactorsMat([1 end],:) = 1;

% 
% import OptMap.KymoAlignment.apply_horizontal_stretching;
spalignedKymo = apply_stretching(alignedkymo, stretchFactorsMat);
% 
% meanFirst = round(nanmean(pathsColIdxs(:,1)));
% for s = 1:size(pathsColIdxs,1)
%     almat(s,:) = circshift(almat(s,:),[0 meanFirst-pathsColIdxs(s,1)]);
%     alignedKymoS(s,:) = circshift(alignedKymoS(s,:),[0 meanFirst-pathsColIdxs(s,1)]);
% end



spalignedBitmask = round(apply_stretching(double(cutMaskF), stretchFactorsMat)) > 0;

if nargout >=5
    f=figure('units','normalized','outerposition',[0 0 1 1])

    tiledlayout(5,1)
    nexttile
    imagesc(kymo)
    title('Original')
    nexttile
    imagesc(alignedkymo)
    title('Shift aligned')

     nexttile
     title('Features using matchpairs')
    hold on
    mat = zeros(size(alignedkymo));
    try
    for s = 1:length(sourcesLong)
        plot( sourcesLong{s}(:, 2),sourcesLong{s}(:, 1),'|','linewidth', 2)
        for k=1:size(sourcesLong{s},1)
        mat(sourcesLong{s}(k, 1), sourcesLong{s}(k, 2)) = 1;
        end
    end
    xlim([1 size(alignedkymo,2)])
    % ylim([3.95, 5.05])
    set(gca, 'YDir','reverse')
    catch
    end
    nexttile
    imagesc(spalignedKymo)
    title('Aligned')
end
% almat = apply_stretching(mat, stretchFactorsMat);



    function [cutKymo,cutMaskF] = shift_align(kymo,bitmask,minOverlap,maxShift)

    % If minOverlap larger than size of kymo, should just return the
    % original kymo.
    
    % shuffle rows
    % shuff = randperm(size(kymo,1));
    shuff = 1:size(kymo,1);
    kymoShuffled = kymo(shuff',:);
    bitmaskShuffled = bitmask;

    lenR = size(kymoShuffled,2);

    import OptMap.KymoAlignment.SPAlign.masked_cc_corr;

    shiftF = zeros(1,size(bitmaskShuffled,1));
    score = zeros(1,size(bitmaskShuffled,1));
    b = 1;

    %% define forward/backward for more accurate 
    for i = 1:size(bitmaskShuffled,1)-1        
        %         import SignalRegistration.masked_cc_corr;
        [ xcorrs, numElts ] = masked_cc_corr([kymoShuffled(i+1,:) nan(1,lenR)] , [kymoShuffled(i,:)  zeros(1,lenR)],[bitmaskShuffled(i+1,:) zeros(1,lenR)],[bitmaskShuffled(i,:)  zeros(1,lenR)],minOverlap ); %todo: include division to k to reduce mem
    %   
        xcorrs = xcorrs./(sqrt(numElts)); % from VALMOD
        
        allowedShifts = b-maxShift:b+maxShift;
        allowedShifts = mod(allowedShifts-1,2*lenR)+1;
        % only keep things within maxShift
        xcorrs2 = nan(size(xcorrs));
        xcorrs2(allowedShifts) = xcorrs(allowedShifts);

        [score(i+1), b] = min(xcorrs2);

        if b> lenR
            b = b-2*lenR;
        end
    %     [b b-2*lenR]
        shiftF(i+1) = shiftF(i) + b -1;
        
    end

    %     addedPxF = max(abs(shiftF));

    %     [ xcorrs, numElts ] = masked_cc_corr([kymoShuffled(1,:) nan(1,lenR)] , [kymoShuffled(end,:)  zeros(1,lenR)],[~isnan(kymoShuffled(1,:)) zeros(1,lenR)],[~isnan( kymoShuffled(end,:))  zeros(1,lenR)],minOverlap ); %todo: include division to k to reduce mem


        alignKymoTemp = [zeros(size(kymoShuffled)) kymoShuffled zeros(size(kymoShuffled))];
        kymoMaskF = [zeros(size(kymoShuffled)) bitmaskShuffled zeros(size(kymoShuffled))];

        alignedKymoF = alignKymoTemp;

        for i=2:size(alignedKymoF,1)
            alignedKymoF(i,:) = circshift(  alignKymoTemp(i,:),[0  shiftF(i)]);
            kymoMaskF(i,:) = circshift(  kymoMaskF(i,:),[0  shiftF(i)]);
        end

        fullColumns = find(sum(kymoMaskF,1)>=1); % in case of outliers this might be inaccurate
        cutKymo = alignedKymoF(:,fullColumns(1):fullColumns(end));
        cutMaskF = kymoMaskF(:,fullColumns(1):fullColumns(end));

        cutMaskF = bwareafilt(logical(cutMaskF),1);


    end

    function [vec,Ws,Ypk,Xpk,Wpk] = find_mat_peaks(alignedkymo, gfiltSigma,maxNumFeaturesSoughtK, minPeakDistance)

        % apply gaussian filt
        filtImg = imgaussfilt(alignedkymo,gfiltSigma);
        % figure,imagesc(filtImg)
        % filtImg = alignedkymo;


        vec =cell(1,size(filtImg,1));
        Ws =cell(1,size(filtImg,1));

        for i=1:size(filtImg,1) % we loose last row
            [Ypk,Xpk,Wpk] =findpeaks(filtImg(i,:),'SortStr','descend','NPeaks',maxNumFeaturesSoughtK,'MinPeakDistance',minPeakDistance); % FW for more accuracy?
            vec{i} = Xpk';
            Ws{i} = Ypk';
        end
    end


    function [sourcesLong,pathsColIdxs,longFeatures,sources] = track_peaks(alignedkymo,vec,Ws,distPar)


        % finding k features
        maxDistFrames = inf;

        % track features
        import OptMap.KymoAlignment.SPAlign.track_sources;
        [sources] = track_sources(vec,distPar,maxDistFrames,Ws);
        % 
        % % consider sources present in half or more frames
        % longFeatures = find(cellfun(@(x) size(x,1),sources)>0.9*size(filtImg,1)-1);
        longFeatures = find(cellfun(@(x) size(x,1),sources)==size(alignedkymo,1));

        sourcesLong = sources(longFeatures);
        % 



        pathsColIdxs = nan(size(alignedkymo,1), length(sourcesLong));
        for i=1:length(sourcesLong)
            pathsColIdxs(sourcesLong{i}(:,1),i) = sourcesLong{i}(:,2);
        end
% % 
%%% If some features are not over all rows, should do something like the
%%% following to fill them up:
% % for j=1:size(pathsColIdxs,2)
% %     % fill up nan's with nearest nonnan
% %     m = flipud(pathsColIdxs(:,j));
% %     t = ~isnan(m);
% %     ii = cumsum(t);
% %     ii(ii==0) = 1;
% %     ii = bsxfun(@plus, [0 ii(end,1:end-1)],ii);
% % 
% %     % for s=1:size(ii,2)
% %     %     ii(:,s) = bsxfun(@plus, [ii(end,s)],ii(:,s));
% %     % end
% %     m1 = m(t);
% %     pathsColIdxs(:,j) = flipud(m1(ii));
% % end

% end



    end
end


function [pathsColIdxs ] = align_old(shiftAlignedKymo,shiftAlignedKymoMask)
    
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
    
    
%     if alignmentSuccessTF
        %The k shortest paths through the (prealigned) kymograph are found
    import OptMap.KymoAlignment.find_k_features;
    pathsColIdxOffsets = find_k_features(smoothImg(:, leftColIdx:rightColIdx), maxFeatureMovementPx, typicalFeatureWidth, maxNumFeaturesSoughtK); 

%     alignmentSuccessTF = not(isempty(pathsColIdxOffsets));
    %         if alignmentSuccessTF
    pathsColIdxOffsets = horzcat(pathsColIdxOffsets{:, 1});
    pathsColIdxs = leftColIdx + pathsColIdxOffsets - 1;

            
%         end
end


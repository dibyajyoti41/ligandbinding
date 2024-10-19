function [bitmask, posY,mat,meanD,varD,badMol] = median_filt_alt(km, filterM,bgSigma,sSigma)

% use median filter to detect edges
%   Args:
%       km - cell of input kymograph with the same intensity settings
%       filterM - filter size for filtering
%
%   Returns:
%       posY - positions left/right
%       mat - output thresholded image
%       meanD  - mean of background
%       varD - variance of background

    if nargin < 2
        filterM = [5 15];
    end
    
    if nargin < 3
        bgSigma = 3;
    end
    
    if nargin < 4 
        sSigma = 5;
    end

    
    vecvals = cell(1,length(km));
    for i=1:length(km)
        vecvals{i} = km{i}(:);
    end
    fullmat = double(cell2mat(vecvals'));
    fullmat(fullmat< 1) = nan; % assume nan's everything that is 0
    fullmat = fullmat(~isnan(fullmat)); % remove nan's
    
    % we can find indices of all elements in specific range to identify
    % outliers. Here we want to find kymo's that have intensities less than
    % 20
%    a= cellfun(@(x) min(x(:)),km);
%    km2 = km;
%    km2(41) = [];
%        vecvals = cell(1,length(km2));
%     for i=1:length(km2)
%         vecvals{i} = km2{i}(:);
%     end
%         fullmat = double(cell2mat(vecvals'));
%     fullmat(fullmat< 1) = nan; % assume nan's everything that is 0
%     fullmat = fullmat(~isnan(fullmat)); % remove nan's

    %     https://se.mathworks.com/help/stats/half-normal-distribution.html
    
    % use half normal distribution
    [dat,vals] = hist(fullmat,unique(fullmat)); % histogram for unique values/ could fit noise dist to the left
    
    
    [pos,pk] = findpeaks(imgaussfilt(dat,sSigma),'NPeaks',1);
    
%     % finding noise peak
%     figure,plot(dat)
%     hold on
%     findpeaks(imgaussfilt(dat,sSigma),'NPeaks',1)
% %     
    leftHalfDist = fullmat(fullmat<=vals(pk));
    
    % params of half norma dist
    % https://en.wikipedia.org/wiki/Truncated_normal_distribution#One_sided_truncation_(of_lower_tail)[5]
    meanD = vals(pk);%mean(leftHalfDist)+std(leftHalfDist)*sqrt(2/pi);
    varD = std(leftHalfDist)^2/(1-2/pi);
    
    numSigmas = 3;
    bg = fullmat(fullmat<meanD+numSigmas*sqrt(varD));
    
% %     vals = kmeans(fullmat,2);
% %     % alternative: thresh based on min point between first two peaks
% % 
% %     % check which is larger
% %     mat1 = fullmat(vals==1);
% %     mat2 = fullmat(vals==2);
% % 
% %     if mean(mat1)<mean(mat2)
% %         bg = mat1;
% %     else
% %         bg = mat2;
% %     end

    
    
    % filter, same number of points as signal filter
    filterS = [filterM(1)*filterM(2) 1];

    % calculate mean and standard deviation
    filteredVals = (medfilt2(double(bg),filterS,'symmetric'));
   [datFilt,valsFilt] = hist(filteredVals,unique(filteredVals)); % histogram for unique values/ could fit noise dist to the left
%     figure,plot(datFilt)
%     hold on
%     findpeaks(imgaussfilt(datFilt,sSigma),'NPeaks',1)
% %
    
    [posFilt,pkFilt] = findpeaks(imgaussfilt(datFilt,sSigma),'NPeaks',1);
    
     leftHalfDistFilt = filteredVals(filteredVals<=valsFilt(pkFilt));

    % these are more or less the same as for unfilt?
    threshval = valsFilt(pkFilt); % this we should also do for left distribution 
    threshstd = sqrt(std(leftHalfDistFilt)^2/(1-2/pi));
    %nanstd(medfilt2(double(bg),filterS,'symmetric'));
    % number of standard deviations to allow
%     bgSigma = 1;

    posY = cell(1,length(km)); % put into function! medfilt based edge detection. Less accurate for single-frame stuff
    mat = cell(1,length(km));
    bitmask =  cell(1,length(km));

% 
% i= 106;
% K = medfilt2(km{i},filterM,'symmetric') > threshval+bgSigma*threshstd;
% figure;   
% imshowpair(imresize(K,[200 500]),imresize(km{i},[200 500]), 'ColorChannels','red-cyan'  )

    
%     
% % %     
% figure,tiledlayout(11,11,'TileSpacing','none','Padding','none')
% for i=1:length(km)
%     i
% %         km(isnan(kymos{1}{i}))=0;
%     K = medfilt2(km{i},filterM,'symmetric') > threshval+bgSigma*threshstd;
%     mat{i} = K;
%     [labeledImage, numBlobs] = bwlabel(K);
%     
%         
%     
% %     figure,tiledlayout(2,1);nexttile
% %     imagesc(K)
% %     nexttile
% %     imagesc(km{i})
%     
%     nexttile 
%     imshowpair(imresize(K,[200 500]),imresize(km{i},[200 500]), 'ColorChannels','red-cyan'  )
%     title(num2str(i));
% end


badMol = zeros(1,length(km));
for i=1:length(km)
    K = medfilt2(km{i},filterM,'symmetric') > threshval+bgSigma*threshstd;
%     K = medfilt2(km{i},filterM,'symmetric') > meanD+bgSigma*varD; % should give correct also for unmedian filtered. Sigma bigger though?

%     figure;    imshowpair(imresize(K,[200 500]),imresize(km{i},[200 500]), 'ColorChannels','red-cyan'  )
%     f

%     figure,imagesc(km{i})
    mat{i} = K;
    [labeledImage, numBlobs] = bwlabel(K);
    props = regionprops(labeledImage, 'Area');
    sortedVals = sort([props.Area],'desc');


    if numBlobs == 1 || (length(sortedVals)>1 && sortedVals(2) < 50)

        % if there are two regions, could split into two molecules.
        [maxArea, largestIndex] = max([props.Area]);
        labK = labeledImage==largestIndex; % either just max or create a loop here
        posY{i}.leftEdgeIdxs = arrayfun(@(x) find(labK(x,:) >0,1,'first'),1:size(labK,1));
        posY{i}.rightEdgeIdxs = arrayfun(@(x) find(labK(x,:) >0,1,'last'),1:size(labK,1)); 
        bitmask{i} = zeros(size(K));

        for j=1:size(bitmask{i},1)
            bitmask{i}(j,posY{i}.leftEdgeIdxs(j):posY{i}.rightEdgeIdxs(j)) = 1;
        end
    else
        posY = [];
        badMol(i)=1;
    end
end

% ix = 2;
% kymo = km{ix};




end


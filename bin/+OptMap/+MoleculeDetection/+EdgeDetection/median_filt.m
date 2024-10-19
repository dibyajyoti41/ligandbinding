function [bitmask, posY,mat] = median_filt(km, filterM,bgSigma)

% use median filter to detect edges
%   Args:
%       km - cell of input kymograph with the same intensity settings
%       filterM - filter size for filtering
%
%   Returns:
%       posY - positions left/right
%       mat - output thresholded image

    if nargin < 2
        filterM = [5 15];
    end
    
    if nargin < 3
        bgSigma = 1;
    end

    vecvals = cell(1,length(km));
    for i=1:length(km)
        vecvals{i} = km{i}(:);
    end
    fullmat = double(cell2mat(vecvals'));
    fullmat(fullmat< 1) = nan; % assume nan's everything that is 0
    fullmat = fullmat(~isnan(fullmat)); % remove nan's

    vals = kmeans(fullmat,2);
    % alternative: thresh based on min point between first two peaks

    % check which is larger
    mat1 = fullmat(vals==1);
    mat2 = fullmat(vals==2);

    if mean(mat1)<mean(mat2)
        bg = mat1;
    else
        bg = mat2;
    end
    
    % filter, same number of points as signal filter
    filterS = [1 filterM(1)*filterM(2)];

    % calculate mean and standard deviation
    threshval = nanmean(medfilt2(bg,filterS,'symmetric'));
    threshstd = nanstd(medfilt2(bg,filterS,'symmetric'));
    % number of standard deviations to allow
%     bgSigma = 1;

posY = cell(1,length(km)); % put into function! medfilt based edge detection. Less accurate for single-frame stuff
mat = cell(1,length(km));
bitmask =  cell(1,length(km));
for i=1:length(km)
%         km(isnan(kymos{1}{i}))=0;
    K = medfilt2(km{i},filterM,'zeros') > threshval+bgSigma*threshstd;
    mat{i} = K;
    [labeledImage, numBlobs] = bwlabel(K);

    props = regionprops(labeledImage, 'Area');
    [maxArea, largestIndex] = max([props.Area]);


    try
        labK = labeledImage==largestIndex; % either just max or create a loop here
        posY{i}.leftEdgeIdxs = arrayfun(@(x) find(labK(x,:) >0,1,'first'),1:size(labK,1));
        posY{i}.rightEdgeIdxs = arrayfun(@(x) find(labK(x,:) >0,1,'last'),1:size(labK,1)); 
        bitmask{i} = zeros(size(K));

        for j=1:size(bitmask{i},1)
            bitmask{i}(j,posY{i}.leftEdgeIdxs(j):posY{i}.rightEdgeIdxs(j)) = 1;
        end
    catch
       posY{i} = [];
    end

end

% ix = 2;
% kymo = km{ix};




end


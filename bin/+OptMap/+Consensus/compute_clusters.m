function [lC, clusterMeanCenters ] = compute_clusters(lens, lenRangeFactor )
    % compute_clusters. Computes clusters for barcode lengths

    %   Args:
    %       barcodes 
    %       consensusSettings
    %
    %   Returns:
    %       lC, clusterMeanCenters

    % written by Albertas Dvirnas
    
    % We want to run the cluster consensus over barcodes of different
    % lengths. If the length is within len/lR and len*lR, we deem to
    % barcodes to belong to the same cluster. 
    
    % range factor
    lR = lenRangeFactor;
    
%     lens = cellfun(@length, barcodes); % barcode lengths
    nB = length(lens); % number barcodes
    % the same, so should be removed.
    
    lC = zeros(nB, 1);
    while any(lC == 0)
        tmp_countInLenRange = zeros(nB, 1);
        tmp_stdForLensInRange = zeros(nB, 1);
        for k = 1:nB
            if lC(k) == 0 % if k'th barcode does not belong to a cluster yet
                lenb = lens(k); % choose current length 
                % compute all barcodes that this would be center of in a
                % cluster
                mask = (lens >= (1/lR).*lenb & (lens <= lR.*lenb)) & (lC == 0); 
                tmp_countInLenRange(k) = sum(mask); % number of barcodes
                tmp_stdForLensInRange(k) = std(lens(mask)); % std
            end
        end
        % take the one that has maximum count
        [~,maxLenIdx] = max(tmp_countInLenRange);  
        % if there are a few indice with the same count, take the one with
        % minimum std
        [~, tmp_mvi] = min(tmp_stdForLensInRange(maxLenIdx));
        k = maxLenIdx(tmp_mvi);

        lenb = lens(k);
        mask = (lens >= (1/lR).*lenb & (lens <= lR.*lenb)) & (lC == 0);
        lC(mask) = max(lC(:)) + 1;
    end
    
    % number of clusters
    tmp_numLenClusters = max(lC);
    
    % average length in each cluster
    clusterMeanCenters = zeros(tmp_numLenClusters, 1);
    for tmp_lenClusterNum = 1:tmp_numLenClusters
        tmp_lenClusterMask = (lC == tmp_lenClusterNum);
        tmp_commonLength_pixels = round(mean(lens(tmp_lenClusterMask)));
        clusterMeanCenters(tmp_lenClusterNum) = tmp_commonLength_pixels;
    end

    

end


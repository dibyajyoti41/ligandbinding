function [ alignedKym,unAlignedKymoMoleculeMask,alignedKymoMoleculeMask,backgroundKym, ssdCoef] = ssd_algn(kymoToAlign,sets)
    % ssd_algn
    % Generates experimental barcodes that are aligned to theory
    %     Args:
    %         kymoToAlign: kymos to be aligned
    %         sets: settings  
    %     Returns:
    %          alignedKym,unAlignedKymoMoleculeMask,alignedKymoMoleculeMask,backgroundKym, ssdCoef

    import OptMap.KymoAlignment.SSDAlign.ssd_algn_compute_values;

    % coefficients to be saved
    ssdCoef.left= [];
    ssdCoef.cor = [];
    ssdCoef.tot = [];
    ssdCoef.shift = [];
    
    % filter size to filter by
    filterSize = sets.bitmasking.psfSigmaWidth_nm/sets.bitmasking.prestretchPixelWidth_nm;
 
    % edge indexes
   leftEdgeIdxs = zeros(1,size(kymoToAlign,1));
   rightEdgeIdxs = zeros(1,size(kymoToAlign,1));

   % initialize kymograph mask
   unAlignedKymoMoleculeMask = zeros(size(kymoToAlign));
   % initialize filtered kymo
   kym = cell(1,size(kymoToAlign,1));
   for i=1:size(kymoToAlign,1)
        % filter the row of a kymo using imgaussfilt and filtersize
        % computed before
        kym{i} = imgaussfilt(kymoToAlign(i,:),filterSize);
        
        % use kmeans to separate indexes in background and forward
        [idx1,~] = kmeans(kym{i}',2);
        % make sure indexes are unique
        [~,~,idx1] = unique(idx1,'stable');
        idx1 = idx1-1;
        %
        % find first signal index
        leftEdgeIdxs(i) = find(idx1,1,'first');
        % find last signal index
        rightEdgeIdxs(i) = find(idx1,1,'last');
        % update unaligned kymograph molecule mask row
        unAlignedKymoMoleculeMask(i,leftEdgeIdxs(i):rightEdgeIdxs(i))= ones(1,length(rightEdgeIdxs(i))-length(leftEdgeIdxs(i))+1);             
   end
    
   % this will store aligned kymograph mask
    alignedKymoMoleculeMask = unAlignedKymoMoleculeMask;

    % first row serves as a references
    kym1mol = kym{1}(leftEdgeIdxs(1):rightEdgeIdxs(1)); 
    % edge pixels
    edgePixels = round(sets.bitmasking.untrustedPx);  
    % bitmask of first molecule, only bitmask the edges. 
    kym1bit = ones(1,length(kym1mol));
    kym1bit(1:edgePixels) = zeros(1,length(edgePixels));
    kym1bit(end-edgePixels+1:end) = zeros(1,length(edgePixels));
    
    % this will store the aligned kymograph
    alignedKym = nan(size(kymoToAlign));
    alignedKym(1,:) = kymoToAlign(1,:);
    
    
    for i=2:size(kymoToAlign,1)
         % the row to compare to
         kym2mol = kym{i}(leftEdgeIdxs(i):rightEdgeIdxs(i)); 
         kym2bit = ones(1,length(kym2mol));
         kym2bit(1:edgePixels) = zeros(1,length(edgePixels));
         kym2bit(end-edgePixels+1:end) = zeros(1,length(edgePixels));
      
         % we compute the ssd values for the two, we limit the shift
         % between the two by "edgepixels"
         [ssdV,~, indices] = ssd_algn_compute_values(kym1mol,kym2mol, kym1bit, kym2bit, edgePixels);
        
         % save the ssd coefficients to ssd struct
         ssdCoef.left = [ ssdCoef.left; ssdV];
         
         % find where the minimum is
         [~,b] = min(ssdV);
         relShift = indices(b);

         % correct place for the comparison of two is store here
         ssdCoef.cor = [ ssdCoef.cor b];
         
         % how much we shift depends on relShift and different between the
         % starting position of the two barcodes. TODO: check if no need
         % for +-1 here
         ssdCoef.shift  =[ssdCoef.shift relShift+leftEdgeIdxs(1)-leftEdgeIdxs(i)];
         % proceed with the shift here
         alignedKym(i,:) = circshift(kymoToAlign(i,:),[0,relShift+leftEdgeIdxs(1)-leftEdgeIdxs(i)]);
         % same for the mask
         alignedKymoMoleculeMask(i,:) = circshift(alignedKymoMoleculeMask(i,:),[0,relShift+leftEdgeIdxs(1)-leftEdgeIdxs(i)]);
    end
    % finally, the background kymograph
    backgroundKym = nan(size(alignedKym));
    backgroundKym(~alignedKymoMoleculeMask) = alignedKym(~alignedKymoMoleculeMask);
end


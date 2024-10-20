function [maxcoef, pos, or, secondPos, lenM,dist] = unmasked_MASS_PCC(longVec, shortVec, shortVecBit,longVecBit,k,isLinearTF,numPixelsAroundBestTheoryMask)
    % masked MASS_PCC where short vector can have bitmask
    %
    %   Args:
    %       shortVec, longVec, shortVecBit,lomgVecBit - short barcode,long barcode, short barcode mask, long barcode bitmask
    %       k - parameter for computing pcc
    %       isLinearTF - whether is linear
    %       numPixelsAroundBestTheoryMask - pixeks around best match to
    %       ignore
    %   Returns:
    %
    %       dist
    
    
    shortVecCut = shortVec(logical(shortVecBit));

    % regular MASS PCC
    [dist] = MASS_PCC([longVec,longVec(1:min(end-1,sum(shortVecBit)-1))],shortVecCut, k);
    
    % shift back
    dist(1,:) = circshift(dist(1,:),[0,1-find(shortVecBit,1,'first')]);
    dist(2,:) = circshift(dist(2,:),[0,1-find(shortVecBit,1,'first')]);
        
    longVecBit = [longVecBit  longVecBit(1:min(end-1,sum(shortVecBit)-1))];
    
    dist(isinf(dist)) = nan; % inf should be undefined.
    import CBT.Hca.UI.Helper.get_best_parameters;
    [maxcoef, pos, or] = get_best_parameters(dist, 3,length(shortVec),isLinearTF,numPixelsAroundBestTheoryMask,longVecBit);
    
    %
    secondPos = find(shortVecBit,1,'first');
    lenM = length(shortVecCut);
    
%     posMax = posMax-secondPos+1;

end


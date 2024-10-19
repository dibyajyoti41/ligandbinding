function [ xcorrs ] = unmasked_pcc_corr( shortVec, longVec, shortVecBit )
    % compute_pcc_nobitmask
	% Computes Pearson correlation coefficient using fft's, no bitmask on
	% the second barcode, the first one just has a simple bitmask
    %     Args:
    %         shortVec: short linear barcode
    %         longVec: long circular barcode
    %         shortVecBit: bitmask of first barcode i.e. 000011110000
    % 
    %     Returns:
    %         xcorrs: PCC values
    %
    shortVecCut = zscore(shortVec(logical(shortVecBit)));

    shortLength = size(shortVecCut,2);
    longLength = size(longVec,2);
    
    shortVecFlip = fliplr(shortVecCut);
     
    longfft = fft(longVec);
    conVec = conj(fft(shortVecCut,longLength));
    % Forward cross correlations
    ccForward = (ifft(longfft.*conVec))/(shortLength-1);

    conVec = conj(fft(shortVecFlip,longLength));
    % Backward cross correlations
    ccBackward = (ifft(longfft.*conVec))/(shortLength-1);
    
    conVec = conj(fft(ones(1,shortLength),longLength));
    movMean =  (ifft(longfft.*conVec))./shortLength;
    movStd =  (ifft(fft(longVec.^2).*conVec));

    stdForward = sqrt((movStd-shortLength.*movMean.^2)./(shortLength-1 ));

    ccForward = ccForward ./ stdForward; 
    ccBackward = ccBackward ./stdForward; % std is the same for both forward and backward case

    xx2 = circshift(ccBackward,[0,1-find(shortVecBit,1,'first')]);
    xx1 = circshift(ccForward,[0,1-find(shortVecBit,1,'first')]);
    xcorrs = [xx1;xx2];
end


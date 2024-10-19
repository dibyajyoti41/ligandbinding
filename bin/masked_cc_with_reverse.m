function [ xcorrs, numElts ] = masked_cc_with_reverse( shortVec, longVec, w1, w2, minLen )
    % masked_cc_with_reverse
    % Computes masked cross correlation using fft's
    %     Args:
    %         shortVec: short linear barcode
    %         longVec: long linear barcode
    %         w1: bitmask of first barcode
    %         w2: bitmask of second barcode
    % 
    %     Returns:
    %         xcorrs: PCC values
    %
    %     Example:
    %       
    % 
    %     % compute X*Y
    % shortLength = 100;
    % longLength = 1000;
    % 
    % shortVec = rand(1,shortLength);
    % longVec = rand(1,longLength);
    %     % bitmasks first bar
    %     w1 = ones(1,shortLength);
    %     w1(1)=0;
    %     w1(end)=0;
    % 
    %     % bitmasks second bar
    %     w2 = ones(1,longLength);
    %     w2(1)=0;
    %     w2(end)=0;
    %   [ xcorrs ] = compute_pcc( shortVec,longVec,w1,w2 )

    if nargin <  5
        minLen = 20;
    end

    longLength = length(longVec);
    shortlen = length(shortVec);
    
    if shortlen > longLength % fix the lengths by adding zeros
        longVec = [longVec zeros(1,shortlen-longLength)];
        w2 = [ w2  zeros(1,shortlen-longLength)];
        longLength = shortlen; 
    end
        
    numElts = zeros(1,longLength);

    % add zeros where bitmask is 0
    shortVec(~w1) = 0;
    longVec(~w2) = 0;

    
    w2fft = fft(w2); % fft
    
    % Compute the number of nonzero elements in each comparison
    conVec = conj(fft(w1,longLength));
    numForward =  round(ifft(w2fft.*conVec));
    numForward(numForward<(minLen)) = nan;
    % compute y^2
    ySq =  (ifft(fft(longVec.^2).*conVec));

    % main part. Compute xiyi
    conVec = conj(fft(shortVec,longLength));
    xiyi = (ifft(fft(longVec).*conVec));

    % now, compute mean of the first vector sq
    conVec = conj(fft(shortVec.^2,longLength));
    xSq =  (ifft(w2fft.*conVec));
    
    ccF = (ySq - 2*xiyi + xSq);

    numElts(1,:) = numForward;

%     % First part, just compute XY
    shortVecFlip = fliplr(shortVec);
    w1Flip = fliplr(w1);     
% 
%     % Compute the number of nonzero elements in each comparison
    conVec = conj(fft(w1Flip,longLength));
    numForward =  round(ifft(w2fft.*conVec));
    numForward(numForward<minLen) = nan;

    ySq =  (ifft(fft(longVec.^2).*conVec));
% 
%     % main part. Compute xiyi
    conVec = conj(fft(shortVecFlip,longLength));
    xiyi = (ifft(fft(longVec).*conVec));

    
%     % now, compute mean of the first vector sq
    conVec = conj(fft(shortVecFlip.^2,longLength));
    xSq =  (ifft(w2fft.*conVec));

    % CCs backward
    ccB = (ySq - 2*xiyi + xSq);

    
    xcorrs = [ccF; ccB];
    numElts(2,:) = numForward;


    

end


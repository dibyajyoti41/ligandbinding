function [ssdV,ssdB, indices] = ssd_algn_compute_values(rawBarcodeFiltered1,rawBarcodeFiltered2, rawBit1, rawBit2, alowedShift)
    % ssd_algn_compute_values
    % Generates ssd values for comparison of two barcodes with allowed
    % shift
    %     Args:
    %         rawBarcodeFiltered1,rawBarcodeFiltered2, rawBit1, rawBit2, alowedShift
    %     Returns:
    %          ssdV,ssdB, indices

    lenDif = length(rawBarcodeFiltered1)-length(rawBarcodeFiltered2);
 
    if lenDif > 0
        rawBarcodeFiltered2 =[rawBarcodeFiltered2 zeros(1,lenDif)];
        rawBit2 = [ rawBit2 zeros(1,lenDif)];
    else
        rawBarcodeFiltered1 =[rawBarcodeFiltered1 zeros(1,-lenDif)];
        rawBit1 = [ rawBit1 zeros(1,-lenDif)];
    end
	ssdV = [];
    ssdB = [];
    indices = -alowedShift:alowedShift ;

    for cS = indices
        indBit = logical(rawBit1.*circshift(rawBit2,[0,cS]));
        r1 = rawBarcodeFiltered1(indBit);
        
        shifted2 = circshift(rawBarcodeFiltered2,[0,cS]);

        r2 = shifted2(indBit);
        ssdV = [ssdV sum((r1-r2).^2/(length(r1)-1))];
    end
    
     rawBarcodeFiltered2 = fliplr(rawBarcodeFiltered2);
     rawBit2 = fliplr(rawBit2);
     
     for cS = indices
        indBit = logical(rawBit1.*circshift(rawBit2,[0,cS]));
        r1 = rawBarcodeFiltered1(indBit);

        shifted2 = circshift(rawBarcodeFiltered2,[0,cS]);
        r2 = shifted2(indBit);
        ssdB = [ssdB sum((r1-r2).^2/(length(r1)-1))];
    end

end


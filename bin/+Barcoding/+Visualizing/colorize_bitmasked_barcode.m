function [barcodeRGB] = colorize_bitmasked_barcode(barcode, barcodeBitmask, inclusionsHueGradients, exclusionsHueGradients, barDims)
    import Barcoding.Visualizing.normalize_intensity;
    import Barcoding.Visualizing.resize_barcode_image;
            
    if not(isequal(size(inclusionsHueGradients), [2, 3])) || any(inclusionsHueGradients(:) < 0) || any(inclusionsHueGradients(:) > 1)
        error('Bad hue');
    end
    if (size(exclusionsHueGradients, 2) ~= 3) || (size(exclusionsHueGradients, 1) ~= 2) || any(exclusionsHueGradients(:) < 0) || any(exclusionsHueGradients(:) > 1)
        error('Bad excluded hue');
    end

    if nargin < 5
        barHeight = 100;
        barWidth = 1;
        barDims = [barHeight, barWidth]; % height and width for each individual bar in the horizontal barcode image
    end
    
    [barcode, barcodeBitmask] = resize_barcode_image(normalize_intensity(barcode, barcodeBitmask), barcodeBitmask, barDims);
    barcodeR = ((exclusionsHueGradients(2,1) - exclusionsHueGradients(1,1))*barcode + exclusionsHueGradients(1,1)).*not(barcodeBitmask) + ((inclusionsHueGradients(2,1) - inclusionsHueGradients(1,1))*barcode + inclusionsHueGradients(1,1)).*barcodeBitmask;
    barcodeG = ((exclusionsHueGradients(2,2) - exclusionsHueGradients(1,2))*barcode + exclusionsHueGradients(1,2)).*not(barcodeBitmask) + ((inclusionsHueGradients(2,2) - inclusionsHueGradients(1,2))*barcode + inclusionsHueGradients(1,2)).*barcodeBitmask;
    barcodeB = ((exclusionsHueGradients(2,3) - exclusionsHueGradients(1,3))*barcode + exclusionsHueGradients(1,3)).*not(barcodeBitmask) + ((inclusionsHueGradients(2,3) - inclusionsHueGradients(1,3))*barcode + inclusionsHueGradients(1,3)).*barcodeBitmask;
    barcodeRGB = cat(3, barcodeR, barcodeG, barcodeB);
end
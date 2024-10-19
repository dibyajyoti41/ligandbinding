function [barcodeImg, barcodeBitmask] = resize_barcode_image(barcodeImg, barcodeBitmask, barDims)
    validateattributes(barcodeImg, {'numeric'}, {'nonnegative', '2d'}, 1);
    validateattributes(barcodeBitmask, {'logical'}, {'nonnegative', '2d'}, 2);
    validateattributes(barDims, {'numeric'}, {'positive', 'integer', 'size', [1 2]}, 3);
    if (not(isequal([1 1], barDims))) % rescale image
        barcodeImg = imresize(barcodeImg, barDims .* size(barcodeImg), 'nearest');
        barcodeBitmask = imresize(barcodeBitmask, barDims .* size(barcodeBitmask), 'nearest');
    end
end
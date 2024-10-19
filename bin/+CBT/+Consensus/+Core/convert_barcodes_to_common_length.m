function [barcodes] = convert_barcodes_to_common_length(rawBarcodes, commonLength)
    %fprintf('Converting barcodes to the same length...\n');
    function [rescaledBarcode] = convert_barcode_to_common_length(rawBarcode)
        v = linspace(1, length(rawBarcode), commonLength);
        rescaledBarcode = interp1(rawBarcode, v);
    end
    barcodes = cellfun(...
        @convert_barcode_to_common_length, ...
        rawBarcodes, ...
        'UniformOutput', false);
end
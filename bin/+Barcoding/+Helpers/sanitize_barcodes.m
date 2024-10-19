function barcodes = sanitize_barcodes(barcodes, bitmasks)
    notCell = false;
    if not(iscell(barcodes))
        notCell = true;
        barcodes = {barcodes};
    end

    if not(iscell(bitmasks))
        notCell = true;
        bitmasks = {bitmasks};
    end

    numBarcodes = length(barcodes);

    numBitmasks = length(bitmasks);
    if not(numBarcodes == numBitmasks)
        error('Number of barcodes and number of bitmasks are inconsistent');
    end

    barcodeLengths = cellfun(@(x) length(x), barcodes);
    barcodeLengths = barcodeLengths(:);
    bitmaskLengths = cellfun(@(x) length(x), bitmasks);
    bitmaskLengths = bitmaskLengths(:);
    if not(isequal(barcodeLengths, bitmaskLengths))
        error('Barcode and bitmask lengths are inconsistent');
    end

    for barcodeNum=1:numBarcodes
        barcodes{barcodeNum}(~bitmasks{barcodeNum}) = NaN;
    end

    if notCell
        barcodes = barcodes{1};
    end
end
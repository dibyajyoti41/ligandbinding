function [consensusInputs] = gen_consensus_inputs_struct(displayNames, rawBarcodes, bpsPerPx_original, clusterScoreThresholdNormalized, commonLength, prestretchUntrustedEdgeLenUnrounded_pixels, prestretchPixelWidth_nm, rawBgs, normSetting)
    rawBarcodeLens = cellfun(@length, rawBarcodes);
    minBarcodePixels = 2;
    if any(rawBarcodeLens < minBarcodePixels)
        error('All barcodes must be at least %d pixels long', minBarcodePixels);
    end
    import CBT.Consensus.Core.convert_barcodes_to_common_length;
    [barcodes] = convert_barcodes_to_common_length(rawBarcodes, commonLength);
    % barcodes = cellfun(@zscore, barcodes, 'UniformOutput', false);

    if strcmp(normSetting, 'background')
        bcInnerFunc = @(bc, bg) (bc - bg);
        barcodeNormalizationFunction = @(x) cellfun(bcInnerFunc, x, rawBgs, 'UniformOutput', 0);
    elseif strcmp(normSetting, 'bgmean')
        bcInnerFunc = @(bc, bg) ((bc - bg) / mean(bc - bg));
        barcodeNormalizationFunction = @(x) cellfun(bcInnerFunc, x, rawBgs, 'UniformOutput', 0);
    else 
        barcodeNormalizationFunction = @(x) cellfun(@zscore, x, 'UniformOutput', 0);
    end

    barcodes = barcodeNormalizationFunction(barcodes);

    stretchFactors = cellfun(@length, rawBarcodes)/commonLength;

    import CBT.Bitmasking.generate_zero_edged_bitmask_row;
    bitmasks = arrayfun( ...
        @(stretchFactor) ...
            generate_zero_edged_bitmask_row(commonLength, round(stretchFactor * prestretchUntrustedEdgeLenUnrounded_pixels)), ...
        stretchFactors(:), ...
        'UniformOutput', false);

    bpsPerPx_stretched = bpsPerPx_original;
    bpsPerPx_stretched(bpsPerPx_stretched > 0) = bpsPerPx_stretched(bpsPerPx_stretched > 0).*bpsPerPx_original(bpsPerPx_stretched > 0);
    numBarcodes = length(barcodes);
    prestretchPixelWidth_nm = prestretchPixelWidth_nm*ones(numBarcodes, 1);
    if size(stretchFactors,2) > size(stretchFactors,1);
        stretchFactors = stretchFactors';
    end
    
    postStretchPixelWidth_nm = prestretchPixelWidth_nm./stretchFactors;
    otherBarcodeData = cell(numBarcodes, 1);
    for barcodeNum=1:numBarcodes
        currBarcodeData = struct;
        currBarcodeData.stretchFactor = stretchFactors(barcodeNum);
        currBarcodeData.bpsPerPx_original = bpsPerPx_original(barcodeNum);
        currBarcodeData.bpsPerPx_stretched = bpsPerPx_stretched(barcodeNum);
        currBarcodeData.nmPerPx_original = prestretchPixelWidth_nm(barcodeNum);
        currBarcodeData.nmPerPx_stretched = postStretchPixelWidth_nm(barcodeNum);
        otherBarcodeData{barcodeNum} = currBarcodeData;
    end
    consensusInputs.otherBarcodeData = otherBarcodeData;
    consensusInputs.barcodes = barcodes;
    consensusInputs.rawBarcodes = rawBarcodes;
    consensusInputs.rawBgs = rawBgs;
    consensusInputs.bitmasks = bitmasks;
    consensusInputs.displayNames = displayNames;
    consensusInputs.clusterScoreThresholdNormalized = clusterScoreThresholdNormalized;
end

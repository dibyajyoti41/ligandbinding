function [consensusStruct, cache] = generate_consensus_for_barcodes(rawBarcodes, displayNames, bpsPerPx_original, barcodeConsensusSettings, cache, rawBgs)
    rawBarcodeLens = cellfun(@length, rawBarcodes);
    minBarcodePixels = 2;
    if any(rawBarcodeLens < minBarcodePixels)
        errMsg = sprintf('All barcodes must be at least %d pixels long', minBarcodePixels);
        error(errMsg);
    end

    if nargin < 2
        displayNames = {};
    end
    if nargin < 3
        bpsPerPx_original = [];
    end
    if nargin < 4
        barcodeConsensusSettings = [];
    end
    if nargin < 5
        cache = [];
    end
    
    
    if isempty(displayNames)
        displayNames = arrayfun(@(barcodeNum) sprintf('barcode %d', barcodeNum), (1:length(rawBarcodes))', 'UniformOutput', false);
    end
    if isempty(bpsPerPx_original)
        import OldDBM.General.Import.prompt_files_bps_per_pixel_wrapper;
        [bpsPerPx_original] = prompt_files_bps_per_pixel_wrapper(displayNames);
    end
    
    consensusStruct = [];
    if isempty(barcodeConsensusSettings)
        rawBarcodeLens = cellfun(@length, rawBarcodes);
        
        import CBT.Consensus.UI.Helper.make_barcode_consensus_settings;
        barcodeConsensusSettings = make_barcode_consensus_settings(rawBarcodeLens);
        if isempty(barcodeConsensusSettings)
            return;
        end
    end
    if isempty(cache) && not(isa(cache, 'containers.Map'))
        cache = containers.Map();
    end
    

    import CBT.Consensus.Import.Helper.gen_consensus_inputs_struct;
    consensusInputs = gen_consensus_inputs_struct(...
        displayNames, ...
        rawBarcodes, ...
        bpsPerPx_original, ...
        barcodeConsensusSettings.clusterScoreThresholdNormalized, ...
        barcodeConsensusSettings.commonLength, ...
        barcodeConsensusSettings.prestretchUntrustedEdgeLenUnrounded_pixels, ...
        barcodeConsensusSettings.prestretchPixelWidth_nm, ...
        rawBgs, ...
        barcodeConsensusSettings.barcodeNormalization ...
    );

    promptToConfirmTF = isfield(barcodeConsensusSettings, 'promptToConfirmTF') && (barcodeConsensusSettings.promptToConfirmTF);
    if promptToConfirmTF
        import CBT.Consensus.Import.confirm_consensus_generation;
        numBarcodes = length(consensusInputs.barcodes);
        aborted = confirm_consensus_generation(numBarcodes);
        if aborted
            fprintf('Aborting consensus input generation\n');
            consensusStruct = [];
            return;
        end
    end

    import CBT.Consensus.Core.make_consensus_as_struct;
    [consensusStruct, cache] = make_consensus_as_struct( ...
        consensusInputs.barcodes, ...
        consensusInputs.bitmasks, ...
        consensusInputs.displayNames,...
        consensusInputs.otherBarcodeData, ...
        consensusInputs.clusterScoreThresholdNormalized, ...
        cache, ...
        consensusInputs.rawBarcodes, ...
        consensusInputs.rawBgs ...
    );
end


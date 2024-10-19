function [consensusStruct, cache] = generate_consensus_for_selected(lm, cache)
    import CBT.Consensus.Import.Helper.generate_barcodes_for_selected_kymos;
    [kymoStructs] = generate_barcodes_for_selected_kymos(lm, true, true);
    if isempty(kymoStructs)
        consensusStruct = [];
        return;
    end

    import CBT.Consensus.Import.Helper.check_kymo_structs_for_consensus_inputs;
    [aborted, displayNames, rawBarcodes, bpsPerPx_original, rawBgs] = check_kymo_structs_for_consensus_inputs(kymoStructs);
    if aborted
        fprintf('Aborting consensus input generation\n');
        consensusStruct = [];
        return;
    end

    rawBarcodeLens = cellfun(@length, rawBarcodes);

    import CBT.Consensus.UI.Helper.make_barcode_consensus_settings;
    barcodeConsensusSettings = make_barcode_consensus_settings(rawBarcodeLens);
    if isempty(barcodeConsensusSettings)
        return;
    end
    barcodeConsensusSettings.promptToConfirmTF = true;
    
    import CBT.Consensus.Core.generate_consensus_for_barcodes;
    [consensusStruct, cache] = generate_consensus_for_barcodes(rawBarcodes, displayNames, bpsPerPx_original, barcodeConsensusSettings, cache, rawBgs);
end

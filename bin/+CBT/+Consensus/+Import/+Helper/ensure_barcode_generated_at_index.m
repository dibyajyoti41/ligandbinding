function [kymoStruct, trueValueList] = ensure_barcode_generated_at_index(kymoIndex, lm, trueValueList, dispSkippingMsg)
    if nargin < 4
        dispSkippingMsg = false;
    end

    kymoDispNames = lm.get_diplay_names(kymoIndex);
    kymoDispName = kymoDispNames{1};

    import CBT.Consensus.Import.Helper.ensure_alignment_at_index;
    [alignedKymo, kymoStruct, trueValueList] = ensure_alignment_at_index(kymoIndex, lm, trueValueList);

    if not(isfield(kymoStruct, 'barcodeGen')) || isempty(kymoStruct.barcodeGen)
        fprintf('Generating barcode for ''%s''...\n', kymoDispName);

        import CBT.Consensus.Import.Helper.gen_barcode_data;
        kymoStruct.barcodeGen = gen_barcode_data(alignedKymo);
    elseif dispSkippingMsg
        fprintf('Skipping barcode generation for ''%s''...\n', kymoDispName);
    end
    trueValueList{kymoIndex} = kymoStruct;
end
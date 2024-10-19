function [alignedKymo, kymoStruct, trueValueList] = ensure_alignment_at_index(kymoIndex, lm, trueValueList, dispSkippingMsg)
    if nargin < 3
        trueValueList = lm.get_true_value_list();
    end
    if nargin < 4
        dispSkippingMsg = false;
    end

    kymoDispName = lm.get_diplay_names(kymoIndex);
    kymoDispName = kymoDispName{1};
    kymoStruct = trueValueList{kymoIndex};
    if not(isfield(kymoStruct, 'alignedKymo')) || isempty(kymoStruct.alignedKymo)
        fprintf('Aligning ''%s''...\n', kymoDispName);
        unalignedKymo = kymoStruct.unalignedKymo;
        import OptMap.KymoAlignment.NRAlign.nralign;
        [alignedKymo] = nralign(unalignedKymo);
        kymoStruct.alignedKymo = alignedKymo;
    else
        alignedKymo = kymoStruct.alignedKymo;
        if dispSkippingMsg
            fprintf('Skipping pre-aligned ''%s''...\n', kymoDispName);
        end
    end
    trueValueList{kymoIndex} = kymoStruct;
end
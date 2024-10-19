function [unfoldedComponents] = unfold_consensus(barcodeKey, barcodeStructsMap)
    import CBT.Consensus.Helper.unfold_consensus;
    barcodeStruct = barcodeStructsMap(barcodeKey);
    if isempty(barcodeStruct.parents)
        unfoldedComponents = {barcodeKey, false, 0};
        return;
    end
    keys = cell(2,1);
    flipTFs = false(2,1);
    circShifts = zeros(2,1);
    unfoldedComponents = cell(2,1);
    for parentComponentNum=1:2
        [keys{parentComponentNum}, flipTFs(parentComponentNum), circShifts(parentComponentNum)] =  barcodeStruct.parents{parentComponentNum}{:};
        unfoldedComponents{parentComponentNum} = unfold_consensus(keys{parentComponentNum}, barcodeStructsMap);
        unfoldedComponents{parentComponentNum}(:,3) =  num2cell([unfoldedComponents{parentComponentNum}{:,3}]' + (1 - 2*[unfoldedComponents{parentComponentNum}{:,2}]')*circShifts(parentComponentNum));
        if flipTFs(parentComponentNum)
            unfoldedComponents{parentComponentNum}(:,2) = cellfun(@not, unfoldedComponents{parentComponentNum}(:,2), 'UniformOutput', false);
        end
    end
    unfoldedComponents = [unfoldedComponents{1}; unfoldedComponents{2}];
end
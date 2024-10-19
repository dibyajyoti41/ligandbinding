function [aborted, displayNames, rawBarcodes, bpsPerPx_original, rawBgs] = check_kymo_structs_for_consensus_inputs(kymoStructs)
    kymoStructs = kymoStructs(:);
    displayNames = cellfun(@(x) x.displayName, kymoStructs, 'UniformOutput', false);
    rawBarcodes = cellfun(@(x) x.barcodeGen.rawBarcode, kymoStructs, 'UniformOutput', false);
    rawBgs = cellfun(@(x) x.barcodeGen.rawBg, kymoStructs, 'UniformOutput', false);
    bpsPerPx_original = cellfun(@(x) x.bpsPerPx, kymoStructs);
    passesFilters = cellfun(@(x) x.passesFilters, kymoStructs);
    filterFailingDisplayNames = displayNames(not(passesFilters));

    aborted = false;
    if (sum(passesFilters) < 2)
        warndlg('Cannot create consensus with less than two barcodes, but there are fewer filter-passing kymos');
        aborted = true;
    end
    if not(any(passesFilters))
        warndlg('All kymographs'' failed to pass filters', 'All Kymos Fail Filters');
        aborted = true;
    end
    if aborted
        displayNames = {};
        rawBarcodes = {};
        bpsPerPx_original = [];
        return;
    end
    if any(not(passesFilters))
        warndlg('Some kymos failed to pass filters and will be excluded from barcode generation', 'Some Kymos Fail Filters');
        disp('Excluded from consensus:');
        disp(filterFailingDisplayNames);
    end
    displayNames = displayNames(passesFilters);
    rawBarcodes = rawBarcodes(passesFilters);
    rawBgs = rawBgs(passesFilters);
    bpsPerPx_original = bpsPerPx_original(passesFilters);
end

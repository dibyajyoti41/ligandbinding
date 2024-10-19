function [notOK, commonLength] = confirm_stretching_is_ok(commonLength, rawBarcodeLens)
    stretchFactors = rawBarcodeLens/commonLength;
    minStretchFactor = min(stretchFactors);
    maxStretchFactor = max(stretchFactors);
    continueGenPrompt = sprintf('Stretch factors to make the barcodes the same length (%d px) will range from %g to %g. Continue?', ...
        commonLength, ...
        minStretchFactor, ...
        maxStretchFactor);
    continueGenChoice = questdlg(continueGenPrompt, 'Continue consensus generation?', 'Yes', 'No', 'Yes');
    notOK = not(strcmp(continueGenChoice, 'Yes'));
end
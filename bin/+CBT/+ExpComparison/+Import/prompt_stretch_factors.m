function [stretchFactors] = prompt_stretch_factors(defaultMaxStretch, defaultStretchIncrement)
    if nargin < 1
        defaultMaxStretch = .1;
    end
    if nargin < 2
        defaultStretchIncrement = .01;
    end

    promptStretchFactorParams = {'Maximum stretch (%):','Stretch increment (%):'};
    defaultStretchParamVals = {num2str(defaultMaxStretch * 100), num2str(defaultStretchIncrement * 100)};
    stretchFactorParamsDlgTitle = 'Stretch Factor Params';
    num_lines = 1;
    options.Resize='on';
    ansStretchParams = inputdlg(promptStretchFactorParams, stretchFactorParamsDlgTitle, num_lines, defaultStretchParamVals,options);
    if isempty(ansStretchParams)
        ansStretchParams = defaultStretchParamVals;
    end
    maxStretch = str2double(ansStretchParams{1})/100;
    stretchIncrement = str2double(ansStretchParams{2})/100;

    if maxStretch > 0
        numStretchFactors = 1 + 2 * floor(maxStretch / stretchIncrement);
        stretchFactors = linspace(1 - maxStretch, 1 + maxStretch, numStretchFactors);
    end
end
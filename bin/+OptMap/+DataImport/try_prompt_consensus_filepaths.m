function [aborted, consensusFilepaths] = try_prompt_consensus_filepaths(promptTitle, multiselectOn)
    if (nargin < 1)
        promptTitle = [];
    end
    if (nargin < 2) || isempty(promptTitle)
        multiselectOn = true;
    end
    if isempty(promptTitle)
        if multiselectOn
            promptTitle = 'Select Consensus Files';
        else
            promptTitle = 'Select Consensus File';
        end
    end
    if multiselectOn
        multiselectStr = 'on';
    else
        multiselectStr = 'off';
    end
    consensusFilepaths = cell(0, 1);

    [consensusFilenames, dirpath] = uigetfile({'*.mat;'; '*.txt;'}, promptTitle, 'MultiSelect', multiselectStr);

    aborted = isequal(dirpath, 0);
    if aborted
        return;
    end
    if not(iscell(consensusFilenames))
        consensusFilenames = {consensusFilenames};
    end
    consensusFilepaths = fullfile(dirpath, consensusFilenames(:));
end
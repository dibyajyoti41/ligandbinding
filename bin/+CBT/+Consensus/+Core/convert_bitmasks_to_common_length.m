function [bitmasks] = convert_bitmasks_to_common_length(rawBitmasks, commonLength)
    %fprintf('Converting bitmasks to the same length...\n');
    function [rescaledBitmask] = convert_bitmask_to_common_length(rawBitmask)
        v = linspace(1, length(rawBitmask), commonLength);
        rescaledBitmask = rawBitmask(round(v));
    end
    bitmasks = cellfun(...
        @convert_bitmask_to_common_length, ...
        rawBitmasks, ...
        'UniformOutput', false);
end



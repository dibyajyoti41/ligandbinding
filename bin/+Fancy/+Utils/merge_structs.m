function [outputStruct] = merge_structs(varargin)
    % MERGE_STRUCTS - merges scalar structs to produce
    %      a new struct with fields of all the structs
    %   each fields value in the outputStruct is the value of the
    %       field for the last struct for which it is defined
    %
    %  inputs:
    %    any number of scalar structs ordered such that fields
    %       values of latter structs are given higher "priority"
    %        such that they override values from earlier structs
    %
    %  output:
    %    a scalar struct containing fields from the input structs
    %       with values for fields defined by the field value for
    %       the last struct to define it
    
    outputStruct = struct;
    cellOfStructs = varargin;
    numStructs = numel(cellOfStructs);
    if numStructs > 0
        for structNum = 1:numStructs
            currStruct = cellOfStructs{structNum};
            validateattributes(currStruct, {'struct'}, {'scalar'}, structNum);
            currFieldnames = fieldnames(currStruct);
            numFields = numel(currFieldnames);
            for fieldNum = 1:numFields
                currFieldname = currFieldnames{fieldNum};
                outputStruct.(currFieldname) = currStruct.(currFieldname);
            end
        end
    end
end
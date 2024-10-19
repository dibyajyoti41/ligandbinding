function [varargout] = extract_fields(aStruct, fieldsToExtract)
    % EXTRACT_FIELDS - extracts the values for the fields specified
    %   from the scalar struct specified
    %
    % Inputs:
    %   aStruct
    %     a struct from which to extract field values from
    %
    %  fieldsToExtract
    %    a cell of strings specifying what fields to extract as
    %    outputs (a single string for a single field may also be
    %    provided)
    %
    % Outputs: (variable)
    %   the values for the fields of the struct in the order they
    %   were provided as inputs to this function
    %
    % Authors:
    %   Saair Quaderi
    
    validateattributes(aStruct, {'struct'}, {'scalar'}, 1);
    
    if ischar(fieldsToExtract)
       fieldsToExtract = {fieldsToExtract};
    elseif not(iscellstr(fieldsToExtract)) || not(isvector(fieldsToExtract))
       error('Expected input to be a cell vector of strings');
    end
    numOutputs = length(fieldsToExtract);
    varargout = cell(1,numOutputs);

    for outputNum = 1:numOutputs
       fieldName = fieldsToExtract{outputNum};
       varargout{outputNum} = aStruct.(fieldName);
    end
end
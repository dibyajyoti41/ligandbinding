function [] = write_tsv(filepath, dataStruct, structFields, columnNames, indices)
    import Fancy.IO.TSV.format_tsv_val;
    
    colEndDelim = '\t';
    rowEndDelim = '\n';
    validateattributes(dataStruct, {'struct'}, {});
    if nargin > 3
        validateattributes(structFields, {'cell'}, {'vector'});
        if not(iscellstr(structFields))
            error('structFields cell array must contain strings');
        end
    else
        structFields = fieldnames(dataStruct);
    end
    if nargin > 4
        validateattributes(columnNames, {'cell'}, {'vector'});
        if not(iscellstr(columnNames))
            error('columnNames cell array must contain strings');
        end
        if length(columnNames) < length(structFields)
            error('columnNames must have the same length as structFields');
        end
    else
        columnNames = structFields;
    end
    if nargin >= 5
        validateattributes(indices, {'numeric'}, {'positive'});
    else
        indices = NaN;
    end
    outCellArray = {};
    if not(all(isfield(dataStruct, structFields)))
        error('struct fields are missing from dataStruct');
    end
    numCols = length(structFields);
    for colNum=1:numCols
        fieldName = structFields{colNum};
        if isempty(fieldName)
            error('empty field name');
        end
        structFieldVal = dataStruct.(fieldName);
        if (not(isempty(structFieldVal)) && not(isvector(structFieldVal)))
            error('struct field value is not a vector');
        end
        if not(iscell(structFieldVal))
            structFieldVal = structFieldVal(:);
            structFieldVal = num2cell(structFieldVal);
        else
            structFieldVal = structFieldVal(:);
        end
        if (colNum == 1)
            if isnan(indices)
                indices = 1:length(structFieldVal);
            end
            maxIndex = max(indices);
            outCellArray = cell(length(indices) + 1, numCols);
            outCellArray(1, :) = columnNames;
        end
        if length(structFieldVal) < maxIndex
            error('Bad index or inconsistent indices');
        end
        outCellArray(2:end, colNum) = structFieldVal(indices); %#ok<AGROW>
    end
    outCellArray = cellfun(@format_tsv_val, outCellArray, 'UniformOutput', false);
    outCellArray = outCellArray.';
    outputFormat = repmat(['%s', colEndDelim], 1, size(outCellArray, 1));
    outputFormat = [outputFormat(1:(end - length(colEndDelim))), rowEndDelim];

    fileID = fopen(filepath, 'w');
    fprintf(fileID, outputFormat, outCellArray{:});
    fclose(fileID);
end

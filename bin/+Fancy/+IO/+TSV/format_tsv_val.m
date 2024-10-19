function strFormattedVal = format_tsv_val(val)
    if ischar(val)
        if contains(val, char(9))
            error('Tabs are disallowed from values');
        end
        strFormattedVal = ['"', strrep(val, '"', '""'), '"'];
        return;
    end
    if isempty(val)
        strFormattedVal = '';
        return;
    end
    if not(isscalar(val))
        error('Value must be either a string or a scalar value');
    end
    if islogical(val)
        if val
            strFormattedVal = 'true';
        else
            strFormattedVal = 'false';
        end
        return;
    end
    if isnumeric(val)
        strFormattedVal = num2str(val);
        return;
    end
    error('Only logical and numeric scalars, empty values, and strings are supported');
end
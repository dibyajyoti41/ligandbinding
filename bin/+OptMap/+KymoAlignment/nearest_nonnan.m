function [filledDataVect, allnanTF] = nearest_nonnan(dataVect, defaultVal)
    if isempty(dataVect)
        validateattributes(dataVect, {'numeric'}, {}, 1);
        filledDataVect = dataVect;
        return;
    else
        validateattributes(dataVect, {'numeric'}, {'vector'}, 1);
    end
    if nargin < 2
        defaultVal = NaN;
    else
        validateattributes(defaultVal, {'numeric'}, {'scalar'}, 2);
    end
    dataVectSz = size(dataVect);
    dataVect = dataVect(:)';
    idxs = 1:length(dataVect);
    nanMask = isnan(dataVect);
    allnanTF = all(nanMask);
    if not(allnanTF)
        if not(any(nanMask))
            filledDataVect = reshape(dataVect, dataVectSz);
            return;
        end
        nonnanMask = not(nanMask);
        nonnanVals = dataVect(nonnanMask);
        if length(nonnanVals) > 1
            filledDataVect = mean([...
                interp1(idxs(nonnanMask), nonnanVals, idxs, 'nearest', 'extrap'); ...
                fliplr(interp1(idxs(fliplr(nonnanMask)), fliplr(nonnanVals), idxs, 'nearest', 'extrap'))], ...
                1); % fill NaNs with the nearest nonnan value (average of nearest nonnan values if nearness is a tie)
            filledDataVect = reshape(filledDataVect, dataVectSz);
            return;
        end

        defaultVal = nonnanVals(1);
    end
    filledDataVect = repmat(defaultVal, dataVectSz);
end
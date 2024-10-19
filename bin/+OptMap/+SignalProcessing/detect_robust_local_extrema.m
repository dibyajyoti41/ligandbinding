function [robustLocalExtremaIdxs, robustLocalExtremaVals] = detect_robust_local_extrema(curve, minValDistBetweenAdjExtrema, enforceNonmaximaExtremaBookends)
    % DETECT_ROBUST_LOCAL_EXTREMA - detects robust local extrema
    %  similar to matlab's findpeaks, but finds both peaks and troughs
    % notes:
    %  far-end pixels can be designated troughs if the pixels
    %    adjacent to them are greater-valued
    %  if pixels adjacent to an extrema are equal, the first pixel
    %    is designated to be the extrema
    %  when enforcing minimum distance, troughs and peaks which 
    %    satisfy the condition are selected greedily leftwards
    %    and rightwards starting from the maximum peak
    %  this method should have an O(n) worst-case time complexity
    %
    % Inputs:
    %   curve
    %     the curve on which to detect extrema
    %   minValDistBetweenAdjExtrema
    %     (optional, defaults to 0)
    %     a minimum value distance between adjacent peaks
    %     and troughs
    %   enforceNonmaximaExtremaBookends
    %     (optional, defaults to true)
    %     true to enforce that the first and last extrema found are
    %     troughs (may return no extrema if this cannot be
    %     satisfied)
    %
    % Outputs:
    %   robustLocalExtremaIdxs
    %     the indices of the robust extrema that are found
    %   robustLocalExtremaVals
    %     the values of the robust extrema that are found
    %
    % Authors:
    %   Saair Quaderi

    validateattributes(curve, {'numeric'}, {'vector', 'real'}, 1);

    if nargin < 2
        minValDistBetweenAdjExtrema = 0;
    else
        validateattributes(minValDistBetweenAdjExtrema, {'numeric'}, {'scalar', 'real', 'nonnegative'}, 2);
    end

    if nargin < 3
        enforceNonmaximaExtremaBookends = true; %first and last extrema are not allowed to be maxima
    else
        validateattributes(enforceNonmaximaExtremaBookends, {'logical'}, {'scalar'}, 3);
    end

    curve = curve(:);
    tmpEnergyLandscape = [NaN; curve; NaN];

    % keep only the first of any adjacent pairs of equal values (including NaN)
    eqAdjValMask = (tmpEnergyLandscape(1:end-1) ~= tmpEnergyLandscape(2:end));
    nonnanMask = ~isnan(tmpEnergyLandscape);
    nonnanMask2 = (nonnanMask(1:end-1) | nonnanMask(2:end));
    eqAdjValMask = eqAdjValMask & nonnanMask2;

    idxs = [1; 1 + find(eqAdjValMask)];

    % take the sign of the first derivs derivative
    derivSigns = sign(diff(tmpEnergyLandscape(idxs)));
    inflectionPointsMask = abs(diff(derivSigns)) > 0;

    % find local extrema
    robustLocalExtremaIdxs = idxs(1 + find(inflectionPointsMask)) - 1;
    robustLocalExtremaVals = curve(robustLocalExtremaIdxs);

    if isempty(robustLocalExtremaIdxs) % no extrema detected (flat)
        return;
    end

    % prevent first and last extrema from being maxima
    %  either by bookending with smaller minima found before first maxima
    %  or after last maxima or, if not possible, by removing the first
    %  maxima, or removing the last maxima from the list of detected
    %  extrema
    if enforceNonmaximaExtremaBookends
        if numel(robustLocalExtremaIdxs) < 2 % assume it is a maxima
            firstIsMaxima = true;
            lastIsMaxima = true;
        else % determine if they are maxima
            firstIsMaxima = robustLocalExtremaVals(1) > robustLocalExtremaVals(2);
            lastIsMaxima = robustLocalExtremaVals(end) > robustLocalExtremaVals(end - 1);
        end
        if firstIsMaxima
            [valFirstMin, idxFirstMin] = min(curve(1:robustLocalExtremaIdxs(1)));
            if valFirstMin == robustLocalExtremaVals(1)
                % remove first maxima
                robustLocalExtremaIdxs(1) = [];
                robustLocalExtremaVals(1) = [];
            else
               % add first minima
               robustLocalExtremaIdxs = [idxFirstMin; robustLocalExtremaIdxs];
               robustLocalExtremaVals = [valFirstMin; robustLocalExtremaVals];
            end
        end
        if lastIsMaxima
            [valLastMin, idxLastMin] = min(curve(robustLocalExtremaIdxs(end):end));
            idxLastMin = robustLocalExtremaIdxs(end) + idxLastMin - 1;
            if valLastMin == robustLocalExtremaVals(end)
                % remove last maxima
                robustLocalExtremaIdxs(end) = [];
                robustLocalExtremaVals(end) = [];
            else
               % add last minima
               robustLocalExtremaIdxs = [robustLocalExtremaIdxs; idxLastMin];
               robustLocalExtremaVals = [robustLocalExtremaVals; valLastMin];
            end
        end
    end

    if minValDistBetweenAdjExtrema > min(abs(diff(robustLocalExtremaVals)))
        % filter out extrema with distances that are too small

        %TODO: the code here could probably be simplified quite a bit

        [tmpMaxVal, tmpMaxIdxIdx] = max(robustLocalExtremaVals);
        extremaInclusion = false(size(robustLocalExtremaIdxs));
        extremaInclusion(tmpMaxIdxIdx) = true;
        leftDone = false;
        tmpLeftMaxVal = tmpMaxVal;
        tmpLeftMaxIdxIdx = tmpMaxIdxIdx;

        rightDone = false;
        tmpRightMaxVal = tmpMaxVal;
        tmpRightMaxIdxIdx = tmpMaxIdxIdx;

        while not(leftDone)
            tmpLeftVals = robustLocalExtremaVals(1:(tmpLeftMaxIdxIdx - 1));
            tmpLeftMins = tmpLeftVals <= (tmpLeftMaxVal - minValDistBetweenAdjExtrema);
            tmpLeftMinIdxIdx = find(tmpLeftMins, 1, 'last');
            if mod(tmpLeftMinIdxIdx, 2) == 0 % a maxima
                tmpLeftMinIdxIdx = tmpLeftMinIdxIdx - 1;
                tmpLeftMinIdxIdx = tmpLeftMinIdxIdx(tmpLeftMinIdxIdx > 0);
            end
            tmpLeftMinVal = robustLocalExtremaVals(tmpLeftMinIdxIdx);
            if isempty(tmpLeftMinIdxIdx)
                leftDone = true;
            else
                extremaInclusion(tmpLeftMinIdxIdx) = true;
                tmpLeftVals = robustLocalExtremaVals(1:(tmpLeftMinIdxIdx - 1));
                tmpLeftMaxs = tmpLeftVals >= (tmpLeftMinVal + minValDistBetweenAdjExtrema);
                tmpLeftMaxIdxIdx = find(tmpLeftMaxs, 1, 'last');
                if mod(tmpLeftMaxIdxIdx, 2) == 1 % a minima
                    tmpLeftMaxIdxIdx = tmpLeftMaxIdxIdx - 1;
                    tmpLeftMaxIdxIdx = tmpLeftMaxIdxIdx(tmpLeftMaxIdxIdx > 0);
                end
                tmpLeftMaxVal = robustLocalExtremaVals(tmpLeftMaxIdxIdx);
                if isempty(tmpLeftMaxIdxIdx)
                    leftDone = true;
                else
                    extremaInclusion(tmpLeftMaxIdxIdx) = true;
                end
            end
        end

        % TODO: the code for the left side could
        %    be made into a function and reused to do this with
        %    just a little cleverness (flipping/index manipulations)
        while(not(rightDone))
            tmpRightVals = robustLocalExtremaVals((tmpRightMaxIdxIdx + 1):end);
            tmpRightMins = tmpRightVals <= (tmpRightMaxVal - minValDistBetweenAdjExtrema);
            [~, tmpRightMinIdxIdx] = find(tmpRightMins, 1, 'first');
            tmpRightMinIdxIdx = tmpRightMinIdxIdx + tmpRightMaxIdxIdx;
            if mod(tmpRightMinIdxIdx, 2) == 0 % a maxima
                tmpRightMinIdxIdx = tmpRightMinIdxIdx + 1;
                tmpRightMinIdxIdx = tmpRightMinIdxIdx(tmpRightMinIdxIdx < numel(robustLocalExtremaVals));
            end
            tmpRightMinVal = robustLocalExtremaVals(tmpRightMinIdxIdx);
            if isempty(tmpRightMinIdxIdx)
                rightDone = true;
            else
                extremaInclusion(tmpRightMinIdxIdx) = true;
                tmpRightVals = robustLocalExtremaVals((tmpRightMinIdxIdx + 1):end);
                tmpRightMaxs = tmpRightVals >= (tmpRightMinVal + minValDistBetweenAdjExtrema);
                [~, tmpRightMaxIdxIdx] = find(tmpRightMaxs, 1, 'first');
                tmpRightMaxIdxIdx = tmpRightMaxIdxIdx + tmpRightMinIdxIdx;
                if mod(tmpRightMaxIdxIdx, 2) == 1 % a minima
                    tmpRightMaxIdxIdx = tmpRightMaxIdxIdx + 1;
                    tmpRightMaxIdxIdx = tmpRightMaxIdxIdx(tmpRightMaxIdxIdx < numel(robustLocalExtremaVals));
                end
                tmpRightMaxVal = robustLocalExtremaVals(tmpRightMaxIdxIdx);
                if isempty(tmpRightMaxIdxIdx)
                    rightDone = true;
                else
                    extremaInclusion(tmpRightMaxIdxIdx) = true;
                end
            end
        end

        extremaInclusion = find(extremaInclusion);
        if enforceNonmaximaExtremaBookends
            minimaIdxs = extremaInclusion(mod(extremaInclusion, 2) == 1);
            firstMinimaIdx = min(minimaIdxs);
            lastMinimaIdx = max(minimaIdxs);
            extremaInclusion = extremaInclusion((extremaInclusion >= firstMinimaIdx) & (extremaInclusion <= lastMinimaIdx));
            if not(isempty(extremaInclusion))
                [~, newFirstMinimaIdx] = min(robustLocalExtremaVals(1:2:firstMinimaIdx));
                newFirstMinimaIdx = newFirstMinimaIdx*2 - 1;
                extremaInclusion(1) = newFirstMinimaIdx;
                [~, newLastMinimaIdx] = min(robustLocalExtremaVals(lastMinimaIdx:2:end));
                newLastMinimaIdx = newLastMinimaIdx*2 - 1;
                newLastMinimaIdx = newLastMinimaIdx + lastMinimaIdx - 1;
                extremaInclusion(end) = newLastMinimaIdx;
            end
        end
    else
        extremaInclusion = 1:numel(robustLocalExtremaIdxs);
    end

    robustLocalExtremaIdxs = robustLocalExtremaIdxs(extremaInclusion);
    robustLocalExtremaVals = robustLocalExtremaVals(extremaInclusion);
end

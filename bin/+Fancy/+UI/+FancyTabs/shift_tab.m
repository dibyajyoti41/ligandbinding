function shift_tab(hTab, shiftAmount, ignoreRangeError)
    % SHIFT_TAB - shifts the location of the tab by some quantity a tab amongst its siblings from its current position
    %
    % Inputs:
    %   hTab:
    %     the handle of the tab
    %   shiftAmount
    %     the amount the tab should be shifted
    %     negative shifts left and positive shifts right
    %     by the absolute value of this value
    %     an infinite magnitude suggests shifting as far as possible, but
    %     otherwise the shiftAmount's magnitude must be within the
    %     possible range unless ignoreRangeError is provided as true
    %   ignoreRangeError(optional...
    %     defaults to true if shiftAmount is Inf and false otherwise)
    %     whether to ignore shift magnitudes which are larger than possible
    %     by treating them as the maximum possible shift in the direction
    %     specified (true) or to throw a validation error for such shift
    %     magnitudes (false)
    %
    % Side-effects:
    %   moves a tab from its current parent to the new parent
    %   and moves the context menu for it to the figure of
    %   the new parent
    %
    % Authors:
    %   Saair Quaderi
    
    import Fancy.UI.FancyTabs.get_tab_shift_range;
    
    
    validateattributes(shiftAmount, {'numeric'}, {'scalar', 'real'}, 2);
    if isfinite(shiftAmount)
        validateattributes(shiftAmount, {'numeric'}, {'integer'}, 2);
    end
    if shiftAmount == 0
        return
    end
    if nargin < 3
        if abs(shiftAmount) == Inf
            ignoreRangeError = true;
        else
            ignoreRangeError = false;
        end
    else
    validateattributes(ignoreRangeError, {'logical'}, {'scalar'}, 3);
        
    end
    
    [maxLeftShiftAmount, maxRightShiftAmount, oldTabIndex, hTabAllSiblings, hParent] = get_tab_shift_range(hTab);
    
    if not(ignoreRangeError)
        validateattributes(shiftAmount, {'numeric'}, {'>=', maxLeftShiftAmount, '<=', maxRightShiftAmount}, 2);
    else
        if shiftAmount < maxLeftShiftAmount
            shiftAmount = maxLeftShiftAmount;
        elseif shiftAmount > maxRightShiftAmount
            shiftAmount = maxRightShiftAmount;
        end
    end
    
    maxIndex = numel(hTabAllSiblings);
    newTabIndex = oldTabIndex + shiftAmount;
    indices = [(1:(oldTabIndex - 1)), ((oldTabIndex + 1):maxIndex)];
    indices = [indices(1:(newTabIndex - 1)), oldTabIndex, indices(newTabIndex:(maxIndex - 1))];
    set(hParent, 'Children', hTabAllSiblings(indices));
end
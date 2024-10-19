function [maxLeftShiftAmount, maxRightShiftAmount, tabIndex, hTabAllSiblings, hParent] = get_tab_shift_range(hTab)
    % GET_TAB_SHIFT_RANGE - gets the maximum left shift and right shift
    %   possible for a tab amongst its siblings from its current position
    %
    % Inputs:
    %  hTab
    %   the handle of the tab
    %
    % Outputs:
    %   maxLeftShiftAmount
    %     the maximum shifting to go left (as a non-positive number)
    %   maxRightShiftAmount
    %     the maximum shifting to go right (as a non-negative number)
    %   tabIndex
    %     the current index of the tab amongst its siblings
    %   hTabAllSiblings
    %     the handles of the parent's children
    %   hParent
    %     the handle of the tab's parent
    %
    % Authors:
    %   Saair Quaderi
    
    import Fancy.UI.FancyTabs.get_tab_family;
    
    [hParent, hTabAllSiblings, tabIndex] = get_tab_family(hTab);
    maxLeftShiftAmount = 1 - tabIndex;
    maxRightShiftAmount = numel(hTabAllSiblings) - tabIndex;
end
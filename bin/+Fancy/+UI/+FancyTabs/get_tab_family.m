function [hTabParent, hTabAllSiblings, tabIndex] = get_tab_family(hTab)
    % GET_TAB_FAMILY - gets the parent of the tab,  all the children of
    %   of the parent (siblings and self), and the index of the tab amongst
    %   the siblings
    %
    % Inputs:
    %   hTab
    %     the handle of the tab
    %
    % Outputs:
    %   hTabParent
    %     the handle of the tab's parent
    %   hTabAllSiblings
    %     the handles of the tab's parent's children
    %   tabIndex
    %     the current index of the tab amongst its siblings
    %
    % Authors:
    %   Saair Quaderi
    
    validateattributes(hTab, {'matlab.ui.container.Tab'}, {'scalar'}, 1);
    
    hTabParent = get(hTab, 'Parent');
    hTabAllSiblings = get(hTabParent, 'Children');
    tabIndex = find(hTabAllSiblings == hTab, 1);
end
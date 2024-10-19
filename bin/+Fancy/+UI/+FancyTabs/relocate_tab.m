function [] = relocate_tab(hTab, hTabNewParent)
    % RELOCATE_TAB - relocates the tab to become a child of the parent
    %   whose handele is provided
    %
    % Inputs:
    %   hTab
    %     a uitab handle
    %   hTabNewParent
    %     the handle for the new parent for the uitab
    %
    % Side-effects:
    %   moves a tab from its current parent to the new parent
    %   and moves the context menu for it to the figure of
    %   the new parent
    %
    % Authors:
    %   Saair Quaderi

    validateattributes(hTab, {'matlab.ui.container.Tab'}, {'scalar'}, 1);
    validateattributes(hTabNewParent, {'matlab.ui.container.TabGroup'}, {'scalar'}, 2);
    
    if iscell(hTabNewParent) && (numel(hTabNewParent) == 1)
        hTabNewParent = hTabNewParent{1};
    end
    hNewFig = ancestor(hTabNewParent, 'figure');
    set(hTab, 'Parent', hTabNewParent);
    hFigContextMenu = get(hTab, 'UIContextMenu');
    set(hFigContextMenu, 'Parent', hNewFig);
end
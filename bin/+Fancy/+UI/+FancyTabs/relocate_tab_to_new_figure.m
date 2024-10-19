function [hNewFig] = relocate_tab_to_new_figure(hTab, figName)
    % RELOCATE_TAB_TO_NEW_FIGURE - relocates the tab to become a child of a
    %   new tabgroup for a new figure
    %
    % Inputs:
    %  hTab
    %    a uitab handle
    %  figName
    %    the name for the new parent figure
    %
    %  hNewFig
    %    the handle for the newly created figure
    %
    % Side-effects:
    %   creates a new figure with a new tabgroup and moves a tab from its
    %   current parent to the new parent and moves the context menu for it
    %   to the figure of the new parent
    %
    % Authors:
    %   Saair Quaderi
    
    validateattributes(hTab, {'matlab.ui.container.Tab'}, {}, 1);
    
    if nargin < 2
        figName = '';
        defaultFigName = get(hTab, 'Title');
        answers = inputdlg({'Enter new figure name:'}, 'New figure', [1 50], {defaultFigName}, 'on');
        if iscell(answers) && not(isempty(answers))
            figName = answers{1};
        end
    end
    
    
    import Fancy.UI.FancyTabs.relocate_tab;
    
    
    hNewFig = figure(...
        'Name', figName, ...
        'Units', 'normalized', ...
        'OuterPosition', [0.05 0.05 0.9 0.9], ...
        'NumberTitle', 'off', ...
        'MenuBar', 'none', ...
        'ToolBar', 'none');
    hTabNewParent = uitabgroup(hNewFig);
    set(hTab, 'Parent', hTabNewParent);
    
    relocate_tab(hTab, hTabNewParent);
end
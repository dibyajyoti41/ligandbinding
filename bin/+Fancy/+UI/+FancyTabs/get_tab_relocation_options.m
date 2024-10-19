function [relocationCallbacks, relocationLabels] = get_tab_relocation_options(hTab)
    % GET_TAB_RELOCATION_OPTIONS - returns functions that, if used, would
    %  relocate the tab that was provided, and some labels for where the
    %  tab would be relocated to (relocation here means changing the tabs
    %  parent)
    %
    % Inputs:
    %   hTab
    %     a uitab handle
    %
    % Outputs:
    %   relocationCallbacks
    %     a cell array of function handles for functions that provide
    %     potential relocations if called (no inputs, no outputs,
    %     just this side effect)
    %   relocationLabels
    %     cell array of strings describing where the function associated 
    %     with the same index in relocationCallbacks will relocate the tab
    %     to if it is used
    %
    % Authors:
    %   Saair Quaderi

    hFigs = findobj('Type', 'figure');
    import Fancy.UI.FancyTabs.get_figure_titles;
    relocationLabels = get_figure_titles(hFigs);

    hFigOld = ancestor(hTab, 'figure');

    tmp = hFigs ~= hFigOld;
    relocationLabels = relocationLabels(tmp);
    hFigs = hFigs(tmp);

    relocationParents = arrayfun(@(hFig) findobj(hFig, '-depth', 1, 'Type', 'uitabgroup'), hFigs, 'UniformOutput', false);

    tmp = cellfun(@(hTabgroups) ~isempty(hTabgroups), relocationParents);
    relocationParents = cellfun(@(relocationParent) relocationParent(1), relocationParents(tmp), 'UniformOutput', false);
    relocationLabels = relocationLabels(tmp);

    import Fancy.Utils.arg_swap;
    callbackify = @arg_swap;
    
    import Fancy.UI.FancyTabs.relocate_tab;
    relocationCallbacks = arrayfun(@(relocationParent) callbackify(@() relocate_tab(hTab, relocationParent{1})), relocationParents, 'UniformOutput', false);

    relocationLabels = [{'New figure'}; relocationLabels];
    import Fancy.UI.FancyTabs.relocate_tab_to_new_figure;
    relocationCallbacks = [{callbackify(@() relocate_tab_to_new_figure(hTab))};  relocationCallbacks];
end
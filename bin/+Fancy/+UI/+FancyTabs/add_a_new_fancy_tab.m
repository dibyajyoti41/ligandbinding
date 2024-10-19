function [hTab, fancyTabStruct] = add_a_new_fancy_tab(hTabgroup, tabTitle)
    % ADD_A_NEW_FANCY_TAB - adds a new fancy tab to a tab group
    %   with a specified title and returns a struct which contains
    %   many handles that may be of use
    %
    % Inputs:
    %   hTabgroup
    %     the handle of the tab group on which to add the
    %     tab as a child
    %   tabTitle
    %     the title for the new tab
    %
    % Outputs:
    %   hTab: the tab itself that was created
    %
    %   fancyTabStruct: a 1x1 struct containing various handles
    %      as values for its fields or as values for the fields in
    %      the substruct in the field shift. The struct contains
    %      the following information:
    %     hTab: the tab itself that was created
    %     hFig: the tabgroup's figure
    %     hFigContextMenu: the context menu on the
    %       figure
    %     hCloseTab: context menu option that allows one to close
    %       the tab
    %     hShiftContextMenu: context menu option that contains
    %       any suboptions to shift the tab left or right
    %     shift.hFarLeft: context menu option that allows one to
    %       shift the tab to the leftmost position in the tabgroup
    %     shift.hLeft: context menu option that allows one to
    %       shift the tab to the left in the tabgroup by one
    %       position
    %     shift.hRight: context menu option that allows one to
    %       shift the tab to the right in the tabgroup by one
    %       position
    %     shift.hFarLeft:  context menu option that allows one to
    %       shift the tab to the rightmost position in the tabgroup
    %     hRelocateContextMenu: context menu  that will contain
    %       any suboptions to relocate the tab
    %    Note that the context menu is setup with a callback so
    %      that options will be enabled/disabled/populated
    %      as appropriate when the context menu is opened
    %      via the callback that is set for it; the relocation
    %      options in particular will be generated dynamically
    %      based on what other tabbed figures are open at the time;
    %      shifting options will be enabled/disabled dynamically
    %      based on how many sibling tabs a tab has on each side
    %
    %  Side-effects:
    %    Adds a tab to the tab group and context menu options associated
    %    with it
    %
    % Authors:
    %   Saair Quaderi

    import Fancy.Utils.arg_swap;
    import Fancy.UI.FancyTabs.close_tab_confirmation;
    import Fancy.UI.FancyTabs.update_tab_context_menu;
    import Fancy.UI.FancyTabs.shift_tab;

    validateattributes(hTabgroup, {'matlab.ui.container.TabGroup'}, {'scalar'}, 1);

    callbackify = @arg_swap;

    tabBackgroundColor = [1, 1, 1];

    fancyTabStruct.hFig = ancestor(hTabgroup, 'figure');
    fancyTabStruct.hFigContextMenu = uicontextmenu(fancyTabStruct.hFig);
    fancyTabStruct.hTab = uitab(...
        'Parent', hTabgroup,...
        'Title', tabTitle,...
        'BackgroundColor', tabBackgroundColor);

    fancyTabStruct.hCloseTab = uimenu(fancyTabStruct.hFigContextMenu,...
        'Label', 'Close Tab',...
        'Callback', callbackify(@() close_tab_confirmation(fancyTabStruct.hTab)));

    fancyTabStruct.hShiftContextMenu = uimenu(...
        'Parent', fancyTabStruct.hFigContextMenu,...
        'Label', 'Shift tab');

    fancyTabStruct.shift.hFarLeft = uimenu(...
        'Parent', fancyTabStruct.hShiftContextMenu,...
        'Label', 'Far left',...
        'Callback', callbackify(@() shift_tab(fancyTabStruct.hTab, -inf)));

    fancyTabStruct.shift.hLeft = uimenu(...
        'Parent', fancyTabStruct.hShiftContextMenu,...
        'Label', 'Left',...
        'Callback', callbackify(@() shift_tab(fancyTabStruct.hTab, -1)));

    fancyTabStruct.shift.hRight = uimenu(...
        'Parent', fancyTabStruct.hShiftContextMenu,...
        'Label', 'Right',...
        'Callback', callbackify(@() shift_tab(fancyTabStruct.hTab, 1)));

    fancyTabStruct.shift.hFarRight = uimenu(...
        'Parent', fancyTabStruct.hShiftContextMenu,...
        'Label', 'Far right',...
        'Callback', callbackify(@() shift_tab(fancyTabStruct.hTab, inf)));

    fancyTabStruct.hRelocateContextMenu = uimenu(...
        'Parent', fancyTabStruct.hFigContextMenu,...
        'Label', 'Relocate tab');

    set(fancyTabStruct.hFigContextMenu, 'Callback', callbackify(@() update_tab_context_menu(fancyTabStruct)));
    set(fancyTabStruct.hTab, 'UIContextMenu', fancyTabStruct.hFigContextMenu);
    hTab = fancyTabStruct.hTab;
end

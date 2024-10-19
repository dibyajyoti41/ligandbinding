classdef TabbedScreen < handle
    % TABBEDSCREEN 
    
    properties (Access = private)
        TabGroupHandle = gobjects(0, 1);
    end
    
    properties
        DataMap = [];
    end
    
    methods
        function [ts] = TabbedScreen(hPanel)
            ts.TabGroupHandle = uitabgroup(hPanel);
        end
        function [hTab] = create_tab(ts, tabTitle, deleteFcn)
            if nargin < 3
                deleteFcn = [];
            end
            import Fancy.UI.FancyTabs.add_a_new_fancy_tab;
            if not(isvalid(ts.TabGroupHandle))
                warning('Tab group was closed, creating new tab figure');
                hFig = figure('Name', tabTitle);
                hPanel = uipanel('Parent', hFig);
                import Fancy.UI.FancyTabs.TabbedScreen;
                ts = TabbedScreen(hPanel);
            end
            
            hTabgroup = ts.TabGroupHandle;
            hTab = add_a_new_fancy_tab(hTabgroup, tabTitle);
            if not(isempty(deleteFcn)) && isa(deleteFcn, 'function_handle')
                set(hTab, 'DeleteFcn', deleteFcn)
            end
        end
        function [hTabs] = create_tabs(ts, tabTitles, deleteFcns)
            if nargin < 3
                deleteFcns = [];
            end
            if isempty(deleteFcns) || isscalar(deleteFcns)
                deleteFcns = repmat({deleteFcns}, size(ts));
            end
            hTabs = cellfun(...
                @(tabTitle, deleteFcn) ...
                    ts.create_tab(tabTitle, deleteFcn), ...
                    tabTitles(:), deleteFcns(:));
        end
        function [hTab] = get_tab_handle(ts, tabNum)
            hTabgroup = ts.TabGroupHandle;
            hChildTabs = findobj(tsDBM.TabGroupHandle, ...
                '-depth', 1, ...
                'Type', 'uitab');
            hChildTabs = setdiff(hChildTabs, hTabgroup);
            hTab = hChildTabs(tabNum);
        end
        function select_tab(ts, tabID)
            if isa(tabID, 'matlab.ui.container.Tab')
                hTab = tabID;
                hTabgroup = get(hTab, 'Parent');
            elseif isnumeric(tabID)
                hTabgroup = ts.TabGroupHandle;
                tabNum = tabID;
                hTab = ts.get_tab_handle(tabNum);
            end
            if isempty(hTabgroup) || not(isvalid(hTabgroup))
                warning('Tabgroup could not be found');
            end
            set(hTabgroup, 'SelectedTab', hTab);
        end
        
        function [hasValidGuiTF] = has_valid_gui(ts)
            hTabgroup = ts.TabGroupHandle;
            hasValidGuiTF = not(isempty(hTabgroup) || not(isvalid(hTabgroup)));
        end
        
        function [] = set_data(ts, key, value)
            if isempty(ts.DataMap)
                ts.DataMap = containers.Map();
            end
            ts.DataMap(key) = value;
        end
        
        function [value, foundValueTF] = get_data(ts, key)
            value = [];
            m = ts.DataMap;
            foundValueTF = false;
            if isa(m, 'containers.Map')
                foundValueTF = isKey(m, key);
                if foundValueTF
                    value = m(key);
                end
            end
        end
    end
    methods (Static)
        function [ts, hFig] = make_tabbed_screen_in_new_fig(figTitle)
            hFig = figure(...
                'Name', figTitle, ...
                'Units', 'normalized', ...
                'OuterPosition', [0 0.05 1 0.95], ...
                'MenuBar', 'none', ...
                'ToolBar', 'none');
            
            hPanel = uipanel('Parent', hFig);
            
            import Fancy.UI.FancyTabs.TabbedScreen;
            ts = TabbedScreen(hPanel);
        end
    end
end
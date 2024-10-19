classdef AppDataPoolViewScreen < handle
    % APPDATAPOOLVIEWSCREEN - Screen to view list of data pools or
    %   list of items in a data pool managed by the instance
    %   of AppDataPoolMgr which is provided upon calling
    %   AppDataPoolMgr.get_instance
    %
    % Authors:
    %   Saair Quaderi
    
    properties (Constant)
        Version = [0, 0, 1];
        
        HideDataPoolPrefix = 'zzz_';
    end
    properties (Access = private)
        TabbedScreen
        
        TabbedScreenNibling % nibling = child of sibling
        
        VerifyBeforeRemovals = true;
        Verbose = false;
        
        DataPoolID = NaN;
    end
    properties (SetAccess = private)
        DataPoolListManager
    end
    properties
        ShowHidden = true;
    end
    
    
    methods
        function [adpvs] = AppDataPoolViewScreen(ts, dataPoolID)
            % APPDATAPOOLVIEWSCREEN - Constructor for APPDATAPOOLVIEWSCREEN
            %   object that creates a screen to view data in the instance
            %   of AppDataPoolMgr which is provided upon calling
            %   AppDataPoolMgr.get_instance
            %
            % Inputs:
            %   ts (optional; defaults to handle for a new tabbed screen
            %      which is created in a new figure with the name
            %      "App Data Pools Viewer")
            %      TabbedScreen in which to create the screen on which
            %      data from AppDataPoolMgr is displayed
            %   dataPoolID (optional; defaults to NaN)
            %     if NaN, a screen is created in TabbedScreen which lists
            %       the data pool IDs
            %     otherwise if a data pool ID, a screen is created in
            %      Tabbed Screen which lists the keys for the items in the
            %      specified data pool
            %  
            % Outputs:
            %   adpvs
            %      the instance of APPDATAPOOLVIEWSCREEN contructed
            %
            % Authors:
            %   Saair Quaderi
            
            if nargin < 1
                figTitle = 'App Data Pools Viewer';
                import Fancy.UI.FancyTabs.TabbedScreen;
                ts = TabbedScreen.make_tabbed_screen_in_new_fig(figTitle);
            end
            if nargin < 2
                dataPoolID = NaN;
            end

            if isnan(dataPoolID)
                tabTitle = 'App Data Pools Viewer';
            else
                tabTitle = sprintf('Data Pool List: %s', dataPoolID);
            end
            
            import Fancy.AppMgr.AppDataPoolViewScreen;
            import Fancy.AppMgr.AppDataPoolMgr;
            
            adpvs.DataPoolID = dataPoolID;
            
            deleteFcn = @(varargin) adpvs.delete();
            [hTabViewer, ~] = AppDataPoolViewScreen.make_list_viewer_panel_in_new_tab(ts, tabTitle, deleteFcn);
            adpvs.TabbedScreen = ts;
            
            import Fancy.UI.FancyList.FancyListMgr;
            dplm = FancyListMgr();
            adpvs.DataPoolListManager = dplm;
            dplm.set_ui_parent(hTabViewer);
            dplm.make_ui_items_listbox();
            
            adpvs.sync_ui_list_data();
            
            import Fancy.UI.FancyList.FancyListMgr;
            import Fancy.UI.FancyList.FancyListMgrBtnSet;
            
            flmbs1 = FancyListMgrBtnSet();
            flmbs1.NUM_BUTTON_COLS = 2;
            flmbs1.add_button(FancyListMgr.make_select_all_button_template());
            flmbs1.add_button(FancyListMgr.make_deselect_all_button_template());
            
            flmbs2 = FancyListMgrBtnSet();
            flmbs2.NUM_BUTTON_COLS = 1;
            flmbs2.add_button(AppDataPoolViewScreen.make_resync_button_template(adpvs));
            
            if isnan(adpvs.DataPoolID)
                flmbs2.add_button(AppDataPoolViewScreen.make_view_pool_entries_button_template(adpvs));
            end
            dplm.add_button_sets(flmbs1, flmbs2);
        end
        
        function delete(adpvs)
            % DELETE - cleans up references to to objects explicitly to
            %   reduce the risk of memory leak when the object is no
            %   longer needed
            %
            % Inputs:
            %   adpvs
            %      the app data pool view screen object
            %    
            % Authors:
            %   Saair Quaderi
            
            if isvalid(adpvs.DataPoolListManager)
                delete(adpvs.DataPoolListManager);
            end
        end
        function [hasDataPoolsMask] = has_data_pool_ids(~, dataPoolIDQueries)
            % HAS_DATA_POOL_IDS - returns a mask where the value for an
            %   index being true signifies that the associated data pool
            %   is presently being managed in the instance of the
            %   AppDataPoolMgr provided upon calling
            %   AppDataPoolMgr.get_instance
            %
            % Inputs:
            %   ~ [ignored]
            %     expected to be the app data pool view screen object
            %   dataPoolIDQueries
            %     list of data pool IDs whose presence are being queried
            %
            % Outputs:
            %   hasDataPoolsMask
            %     mask where the value for an index being true signifies
            %     that the associated data pool is present, and false
            %     signified that it is not present
            %    
            % Authors:
            %   Saair Quaderi
            
            import Fancy.AppMgr.AppDataPoolMgr;
            [dataPoolIDs] = AppDataPoolMgr.get_pool_ids();
            n = length(dataPoolIDQueries);
            hasDataPoolsMask = false(n, 1);
            [~, idxQueries, ~] = intersect(dataPoolIDQueries, dataPoolIDs);
            hasDataPoolsMask(idxQueries) = true;
        end
        
        function [] = sync_ui_list_data(adpvs)
            % SYNC_UI_LIST_DATA - syncs the displayed list of items with
            %   the state of the data being managed by the instance of the
            %   AppDataPoolMgr provided upon calling
            %   AppDataPoolMgr.get_instance
            %
            % Inputs:
            %   adpvs
            %      the app data pool view screen object
            %    
            % Authors:
            %   Saair Quaderi
            
            import Fancy.AppMgr.AppDataPoolViewScreen;
            dplm = adpvs.DataPoolListManager;
            [itemDisplayNames, dataPoolIDs] = adpvs.get_pool_data(adpvs.DataPoolID);

            numItems = length(dataPoolIDs);
            nonhiddenDataPoolIDIdxs = (1:numItems)';
            if not(adpvs.ShowHidden)
                hideDataPoolPrefix = AppDataPoolViewScreen.HideDataPoolPrefix;
                hideDataMask = strncmp(hideDataPoolPrefix, dataPoolIDs, length(hideDataPoolPrefix));
                nonhiddenDataPoolIDIdxs = nonhiddenDataPoolIDIdxs(not(hideDataMask));
            end
                
            dplm.set_list_items(itemDisplayNames(nonhiddenDataPoolIDIdxs), dataPoolIDs(nonhiddenDataPoolIDIdxs)); 
        end
    end
    
    methods(Static, Access = private)
        function [itemDisplayNamesForListMgr, itemValuesForListMgr] = get_pool_data(dataPoolID)
            % GET_POOL_DATA - retrieves the list of display names and
            %   some associated values for the list
            %
            % Inputs:
            %   dataPoolID
            %     if NaN, information about the list of data pool
            %       is returned (display names include each pool id and the
            %       number of entries in the pool; values are just the pool
            %       IDs)
            %     otherwise if a data pool ID, information about the list
            %        of items in the specified data pool are returned
            %        (display names are the data item ids; values are the
            %        data items associated with the ids)
            %  
            % Outputs:
            %   itemDisplayNamesForListMgr
            %     display names to be listed
            %   itemValuesForListMgr
            %     values associated with the display names
            %
            % Authors:
            %   Saair Quaderi
            
            import Fancy.AppMgr.AppDataPoolMgr;
            if (nargin < 1) 
                dataPoolID = NaN;
            end 
            if isnan(dataPoolID)
                [dataPoolIDs, dataPoolsItemCounts] = AppDataPoolMgr.get_pool_ids();
                numItems = length(dataPoolIDs);
                
                itemValuesForListMgr = dataPoolIDs;
                itemDisplayNamesForListMgr = arrayfun(...
                    @(dataPoolNum) ...
                        sprintf('%s [%d Entries]', ...
                            dataPoolIDs{dataPoolNum}, ...
                            dataPoolsItemCounts(dataPoolNum)), ...
                    (1:numItems)', ...
                    'UniformOutput', false);
                return;
            end
            [dataItems, dataItemIDs] = AppDataPoolMgr.get_data_items(dataPoolID);
            if not(isempty(dataItemIDs)) && not(ischar(dataItemIDs{1}))
                dataItemIDs = cellfun(...
                    @num2str, ...
                    dataItemIDs, ...
                    'UniformOutput', false);
            end
            itemDisplayNamesForListMgr = dataItemIDs;
            itemValuesForListMgr = dataItems;
        end
        
        function [btnLoadMovies] = make_resync_button_template(adpvs)
            % MAKE_RESYNC_BUTTON_TEMPLATE - generates the FancyListMgrBtn
            %   button template that would trigger resyncing the listed
            %   items with the data
            %
            % Inputs:
            %   adpvs
            %      the app data pool view screen object
            %  
            % Outputs:
            %   btnLoadMovies
            %     the FancyListMgrBtn object that serves as a button
            %     template and can be added to FancyListMgrBtnSet which
            %     can be added to the fancy list interface in FancyListMgr
            %     and triggers the list data resyncing callback
            %
            % Authors:
            %   Saair Quaderi
            
            import Fancy.UI.FancyList.FancyListMgrBtn;
            btnLoadMovies = FancyListMgrBtn('Refresh  List', @(~, ~, flm) adpvs.sync_ui_list_data());
        end
        
        function [btnViewPoolEntries] = make_view_pool_entries_button_template(adpvs)
            % MAKE_VIEW_POOL_ENTRIES_BUTTON_TEMPLATE - generates the
            %   FancyListMgrBtn button template that would trigger
            %   displaying the view screen for a selected data pool
            %
            % Inputs:
            %   adpvs
            %      the app data pool view screen object
            %  
            % Outputs:
            %   btnViewPoolEntries
            %     the FancyListMgrBtn object that serves as a button
            %     template and can be added to FancyListMgrBtnSet which
            %     can be added to the fancy list interface in FancyListMgr
            %     and triggers the addition of a data pool view screen for
            %     the data pool entry in the current list manager that has
            %     been selected 
            %
            % Authors:
            %   Saair Quaderi
            
            import Fancy.UI.FancyList.FancyListMgrBtn;
            function [] = try_view_pool_entries_lists_in_tabs(adpvs)
                import Fancy.AppMgr.AppDataPoolViewScreen;
                flm = adpvs.DataPoolListManager;
                selectedItems = flm.get_selected_list_items();
                selectedDataPoolIDs = selectedItems(:, 2);
                numDataPools = size(selectedDataPoolIDs);
                
                if numDataPools == 0
                    return;
                end
                ts2 = adpvs.TabbedScreenNibling;
                if isempty(ts2)
                    ts = adpvs.TabbedScreen;
                    hTabViewer2 = ts.create_tab('Pool Entry Lists');
                    hPanelViewer2 = uipanel('Parent', hTabViewer2);
                    import Fancy.UI.FancyTabs.TabbedScreen;
                    ts2 = TabbedScreen(hPanelViewer2);
                    adpvs.TabbedScreenNibling = ts2;
                end
                for dataPoolNum = 1:numDataPools
                     dataPoolID = selectedDataPoolIDs{dataPoolNum};
                     AppDataPoolViewScreen(ts2, dataPoolID)
                end
                
            end
            btnViewPoolEntries = FancyListMgrBtn('View Pool Entry Lists', @(~, ~, flm) try_view_pool_entries_lists_in_tabs(adpvs));
        end
        
        function [hTabViewer, hPanelViewer] = make_list_viewer_panel_in_new_tab(ts, tabTitle, deleteFcn)
            % MAKE_VIEW_POOL_ENTRIES_BUTTON_TEMPLATE - generates the
            %   FancyListMgrBtn button template that would trigger
            %   displaying the view screen for a selected data pool
            %
            % Inputs:
            %   ts
            %     the Tabbed Screen for which to create the tab
            %   tabTitle
            %     the title for the tab to be created
            %   deleteFcn
            %     the callback function to run when the created tab is
            %     deleted (i.e. closed)
            %  
            % Outputs:
            %   hTabViewer
            %     the handle for the created tab
            %   hPanelViewer
            %     the handle for the panel in the created tab
            %
            % Authors:
            %   Saair Quaderi
            
            hPanelViewer = ts.create_tab(tabTitle, deleteFcn);
            hTabViewer = uipanel('Parent', hPanelViewer);
        end
    end
end 
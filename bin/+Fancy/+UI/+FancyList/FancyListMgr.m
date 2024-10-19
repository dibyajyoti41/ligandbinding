classdef FancyListMgr < handle

    % FANCYLISTMGR - Fancy List Manager
    properties (Constant)
        Version = [0, 0, 1];
        UIStatePoolIDPrefix = sprintf('%suiState_FLM_', Fancy.AppMgr.AppDataPoolViewScreen.HideDataPoolPrefix);
        ListDispNamesFromIdxDataPoolIDSuffix = '_DispNamesFromIdx';
        ListValsFromIdxDataPoolIDSuffix = '_ValsFromIdx';
    end
    properties
        ParentHandle = gobjects(0)
        ItemsListbox = gobjects(0)
        
        
        ButtonSets = cell(0, 1);
        ConfigSettings = Fancy.UI.FancyList.FancyListMgrConfig();
    end
    properties (SetAccess = private)
        UUID
        ListDispNamesFromIdxDataPoolID
        ListValsFromIdxDataPoolID
    end
    
    methods
        function [flm] = FancyListMgr(hParent, configSettings)
            import Fancy.AppMgr.AppDataPoolViewScreen;
            import Fancy.UI.FancyList.FancyListMgrConfig;
            import Fancy.UI.FancyList.FancyListMgr;
            if nargin >= 1
                iptaddcallback(hParent, 'SizeChangedFcn', @(~, ~) flm.reposition_controls());
                flm.ParentHandle = hParent;
            end
            if nargin >= 2
                flm.ConfigSettings = configSettings;
            end
            flm.UUID = char(java.util.UUID.randomUUID);
            flm.ListDispNamesFromIdxDataPoolID = [flm.UIStatePoolIDPrefix, flm.UUID, flm.ListDispNamesFromIdxDataPoolIDSuffix];
            flm.ListValsFromIdxDataPoolID = [flm.UIStatePoolIDPrefix, flm.UUID, flm.ListValsFromIdxDataPoolIDSuffix];
            
            import Fancy.AppMgr.AppDataPoolMgr;
            AppDataPoolMgr.create_new_data_pool(flm.ListDispNamesFromIdxDataPoolID, 'uint64', 'any');
            AppDataPoolMgr.create_new_data_pool(flm.ListValsFromIdxDataPoolID, 'uint64', 'any');
            
        end

        function delete(flm) % Clean up data super explicitly to reduce any memory leak risk
            flm.remove_list_items(flm.get_all_indices());
            import Fancy.AppMgr.AppDataPoolMgr;
            AppDataPoolMgr.clear_data_pool(flm.ListDispNamesFromIdxDataPoolID);
            AppDataPoolMgr.clear_data_pool(flm.ListValsFromIdxDataPoolID);
            AppDataPoolMgr.remove_data_pool(flm.ListDispNamesFromIdxDataPoolID);
            AppDataPoolMgr.remove_data_pool(flm.ListValsFromIdxDataPoolID);
        end
        
        function [] = set_ui_parent(flm, hParent)
            hParentOld = flm.ParentHandle;
            if isa(hParentOld, 'matlab.graphics.GraphicsPlaceholder') || not(isvalid(hParentOld))
                iptaddcallback(hParent, 'SizeChangedFcn', @(~, ~) flm.reposition_controls());
                flm.ParentHandle = hParent;
            else
                error('List manager already has a valid parent handle');
            end
        end
        
        function [] = make_ui_parent(flm)
            hParent = flm.ParentHandle;
            if isa(hParent, 'matlab.graphics.GraphicsPlaceholder')
                hFig = figure('Name', 'List Manager');
                hPanel = uipanel('Parent', hFig);
                hParent = hPanel;
                flm.set_ui_parent(hParent)
                flm.ParentHandle = hParent;
            end
        end
        
        function [] = make_ui_items_listbox(flm)
            hParent = flm.ParentHandle;
            if isa(hParent, 'matlab.graphics.GraphicsPlaceholder')
                flm.make_ui_parent();
                hParent = flm.ParentHandle;
            end
            
            itemsListboxHandle = flm.ItemsListbox;
            if isa(itemsListboxHandle, 'matlab.graphics.GraphicsPlaceholder')
                itemsListboxHandle = uicontrol('Parent', hParent, 'Style', 'listbox', 'Max', Inf, 'Min', 0);
                flm.ItemsListbox = itemsListboxHandle;
            end
            flm.reposition_listbox_control();
        end
        
        function [] = make_ui_buttons(flm)
            hParent = flm.ParentHandle;
            if isa(hParent, 'matlab.graphics.GraphicsPlaceholder')
                flm.make_ui_parent();
            end
            
            numButtonSets = length(flm.ButtonSets);
            for buttonSetNum=1:numButtonSets
                flmbs = flm.ButtonSets{buttonSetNum};
                numButtonsInSet = length(flmbs.Buttons);
                for buttonSetButtonNum=1:numButtonsInSet
                    flmb = flmbs.Buttons{buttonSetButtonNum};
                    flmb.instantiate(flm);
                end
            end
            flm.reposition_buttons();
        end
        
        function [] = reposition_listbox_control(flm)
            import Fancy.UI.FancyPositioning.set_at_pos_nrm_in_px;
            
            config = flm.ConfigSettings;
            
            LEFT_STARTING_POSITION_PX = config.LEFT_STARTING_POSITION_PX;
            TOP_STARTING_POSITION_PX = config.TOP_STARTING_POSITION_PX;
            DEFAULT_PADDING_PX = config.DEFAULT_PADDING_PX;
            MIN_LISTBOX_HEIGHT_PX = config.MIN_LISTBOX_HEIGHT_PX;
            LISTBOX_WIDTH_PX = config.LISTBOX_WIDTH_PX;
            
            % TODO: simplify
            set_at_pos_nrm_in_px(flm.ItemsListbox, [0, 1, 1, 0.25],  ...
                @(posPx) [ ...
                    posPx(1) + LEFT_STARTING_POSITION_PX + DEFAULT_PADDING_PX, ...
                    posPx(2) - max(posPx(4), MIN_LISTBOX_HEIGHT_PX) - DEFAULT_PADDING_PX - TOP_STARTING_POSITION_PX, ...
                    LISTBOX_WIDTH_PX, ...
                    max(posPx(4), MIN_LISTBOX_HEIGHT_PX) ...
                ]);
        end
        
        function [] = reposition_buttons(flm)
            import Fancy.UI.FancyPositioning.set_at_pos_nrm_in_px;
            
            config = flm.ConfigSettings;
            LEFT_STARTING_POSITION_PX = config.LEFT_STARTING_POSITION_PX;
            TOP_STARTING_POSITION_PX = config.TOP_STARTING_POSITION_PX;
            DEFAULT_PADDING_PX = config.DEFAULT_PADDING_PX;
            MIN_LISTBOX_HEIGHT_PX = config.MIN_LISTBOX_HEIGHT_PX;
            
            currButtonRowNum = 1;
            numButtonSets = length(flm.ButtonSets);
            for buttonSetNum=1:numButtonSets
                flmbs = flm.ButtonSets{buttonSetNum};
                
                BUTTON_HEIGHT_PX = flmbs.BUTTON_HEIGHT_PX;
                BUTTON_WIDTH_PX = flmbs.BUTTON_WIDTH_PX;
                NUM_BUTTON_COLS = flmbs.NUM_BUTTON_COLS;
                allSetButtonHandles = cellfun(@(x) x.Handle, flmbs.Buttons, 'UniformOutput', false);
                allSetButtonHandles = reshape(allSetButtonHandles, fliplr([ceil(numel(allSetButtonHandles)/NUM_BUTTON_COLS), NUM_BUTTON_COLS]))';
                
                for currSetRow = 1:size(allSetButtonHandles, 1)
                    buttonTopDistFromTopBottonPx = currButtonRowNum*BUTTON_HEIGHT_PX + (currButtonRowNum + 1)*DEFAULT_PADDING_PX + TOP_STARTING_POSITION_PX;
                    for currSetCol = 1:size(allSetButtonHandles, 2)
                        hButton = allSetButtonHandles{currSetRow, currSetCol};
                        if not(isa(hButton, 'matlab.graphics.GraphicsPlaceholder'))
                            buttonLeftPx = LEFT_STARTING_POSITION_PX + (currSetCol - 1)*BUTTON_WIDTH_PX + currSetCol*DEFAULT_PADDING_PX;
                            % TODO: simplify
                            set_at_pos_nrm_in_px(hButton, [0, 1, 1, 0.25], ...
                                @(posPx) [ ...
                                    buttonLeftPx, ...
                                    posPx(2) - (max(posPx(4), MIN_LISTBOX_HEIGHT_PX) + buttonTopDistFromTopBottonPx), ...
                                    BUTTON_WIDTH_PX, ...
                                    BUTTON_HEIGHT_PX ...
                                ]);
                        end
                    end
                    currButtonRowNum = currButtonRowNum + 1;
                end
            end
        end
        
        function [] = reposition_controls(flm)
            flm.reposition_listbox_control();
            flm.reposition_buttons();
        end
        
        
        function [] = add_button_sets(flm, varargin)
            buttonSets = varargin;
            
            buttonSets = cellfun(@(buttonSet) flm.fix_button_set_props(buttonSet), buttonSets, 'UniformOutput', false);
            flm.ButtonSets = [flm.ButtonSets; buttonSets(:)];
            flm.make_ui_buttons();
        end
        
        function [buttonSet] = fix_button_set_props(flm, buttonSet)
            configSettings = flm.ConfigSettings;
            eachGapPaddingWidthPx = configSettings.DEFAULT_PADDING_PX;
            totalWidthPx = configSettings.LISTBOX_WIDTH_PX;
            numButtonCols = buttonSet.NUM_BUTTON_COLS;
            totalPaddingWidthPx = (numButtonCols - 1)*eachGapPaddingWidthPx;
            maxTotalButtonWidthPx = totalWidthPx - totalPaddingWidthPx;
            individualButtonWidthPx = floor(maxTotalButtonWidthPx/numButtonCols);
            buttonSet.BUTTON_WIDTH_PX = individualButtonWidthPx;
        end
        
        function [trueValueList] = get_true_value_list(flm)
            import Fancy.AppMgr.AppDataPoolMgr;
            [trueValueList, dataItemIDs] =  AppDataPoolMgr.get_data_items(flm.ListValsFromIdxDataPoolID);
            trueValueList = trueValueList(:);
            itemIdxs = cellfun(@double, dataItemIDs(:));
            
            if any(diff(itemIdxs) ~= 1)
                [itemIdxs, sortOrder] = sort(itemIdxs);
                trueValueList = trueValueList(sortOrder(:));
            end
            
            if not(isequal(itemIdxs, flm.get_all_indices())) 
                error('Missing items');
            end
        end
        
        function [displayNames] = get_diplay_names(flm, indices)
            displayNames = get(flm.ItemsListbox, 'String');
            if nargin > 1
                displayNames = displayNames(indices);
            end
            displayNames = displayNames(:);
        end
        

        function [displayValuePairList, displayNames, trueValueList] = get_display_value_pair_list(flm, indices)
            if nargin < 2
                indices = flm.get_all_indices();
            end
            displayNames = flm.get_diplay_names();
            if isempty(displayNames)
                displayValuePairList = cell(0, 2);
                return;
            end
            trueValueList = flm.get_true_value_list();
            displayValuePairList = [displayNames(indices), trueValueList(indices)];
        end
        
        
        function [selectedIndices] = get_selected_indices(flm)
            selectedIndices = sort(get(flm.ItemsListbox, 'Value'), 'descend');
        end
        
        function [allIndices] = get_all_indices(flm)
            if not(isvalid(flm.ItemsListbox))
                allIndices = [];
                return;
            end
            allIndices = (1:length(get(flm.ItemsListbox, 'string')))';
        end
        
        
        function [selectedItems, selectedIndices] = get_selected_list_items(flm)
            selectedIndices = flm.get_selected_indices();
            selectedItems = flm.get_display_value_pair_list(selectedIndices);
        end
        
        function [items] = get_all_list_items(flm)
            items = flm.get_display_value_pair_list(flm.get_all_indices());
        end


        function [] = select_some(flm, newSelectionIndices)
            currentSelectionValues = get(flm.ItemsListbox, 'Value');
            newSelectionIndices = sort(unique([newSelectionIndices(:); currentSelectionValues(:)]), 'descend');
            set(flm.ItemsListbox, 'Value', newSelectionIndices');
        end
        
        function [] = select_all(flm)
            set(flm.ItemsListbox, 'Value', 1:length(get(flm.ItemsListbox, 'String')));
        end
        
        function [] = deselect_all(flm)
            set(flm.ItemsListbox, 'Value', []);
        end
        

        function [] = set_list_items(flm, textItems, valueItems)
            flm.deselect_all();
            textItems = textItems(:);
            valueItems = valueItems(:);
            numItems = length(textItems);
            
            import Fancy.AppMgr.AppDataPoolMgr;
            dataItemIDs = arrayfun(@uint64, (1:numItems)', 'UniformOutput', false);
            AppDataPoolMgr.clear_data_pool(flm.ListDispNamesFromIdxDataPoolID);
            AppDataPoolMgr.clear_data_pool(flm.ListValsFromIdxDataPoolID);
            AppDataPoolMgr.update_data_items(flm.ListDispNamesFromIdxDataPoolID, dataItemIDs, textItems)
            AppDataPoolMgr.update_data_items(flm.ListValsFromIdxDataPoolID, dataItemIDs, valueItems)
            
            set(flm.ItemsListbox, 'String', textItems);
        end
        
        function [] = add_list_items(flm, newTextItems, newValueItems)
            newTextItems = newTextItems(:);
            newValueItems = newValueItems(:);
            numNewListItems = length(newTextItems);
            if (length(newValueItems) ~= numNewListItems)
                error('Length mismatch between list item texts and values');
            end

            oldTextItems = get(flm.ItemsListbox, 'String');
            oldValueItems = flm.get_true_value_list();
            flm.set_list_items([oldTextItems(:); newTextItems(:)], [oldValueItems; newValueItems(:)]);
        end
        
        function [] = add_list_item(flm, newTextItem, newValueItem)
            flm.add_list_items({newTextItem}, {newValueItem});
        end

        function [] = remove_list_items(flm, indicesForRemoval)
            if not(isvalid(flm.ItemsListbox))
                return;
            end
            
            strings = get(flm.ItemsListbox, 'String');
            maxItemIndex = length(strings);
            indicesForRemoval = indicesForRemoval(indicesForRemoval <= maxItemIndex);

            
            remainingIndices = setdiff((1:maxItemIndex)', indicesForRemoval(:), 'stable');
            
            oldTextItems = get(flm.ItemsListbox, 'String');
            oldValueItems = flm.get_true_value_list();
            flm.set_list_items(oldTextItems(remainingIndices), oldValueItems(remainingIndices));
        end
        
        function [] = remove_selected_items(flm)
            selectedIndicesForRemoval = flm.get_selected_indices();
            flm.remove_list_items(selectedIndicesForRemoval);
        end
    end
    
    methods(Static)
        function [flmbsDefault] = make_default_buttons_set(fn_prompt_for_items, fn_on_items_removal, warnBeforeRemovals)
            if nargin < 1
                fn_prompt_for_items = @nonprompt_for_items;
                disableAddItemsButton = true;
            else
                disableAddItemsButton = false;
            end
            if nargin < 2
                fn_on_items_removal = @nop;
            end
            if nargin < 3
                warnBeforeRemovals = true;
            end
            
            import Fancy.UI.FancyList.FancyListMgr;
            import Fancy.UI.FancyList.FancyListMgrBtnSet;
            
            flmbsDefault = FancyListMgrBtnSet();
            flmbsDefault.NUM_BUTTON_COLS = 2;
            flmbsDefault.add_button(FancyListMgr.make_select_all_button_template());
            flmbsDefault.add_button(FancyListMgr.make_deselect_all_button_template());
            flmbsDefault.add_button(FancyListMgr.make_add_items_button_template(fn_prompt_for_items, disableAddItemsButton));
            flmbsDefault.add_button(FancyListMgr.make_try_remove_items_button_template(warnBeforeRemovals, fn_on_items_removal));
            
            function [aborted, displayNames, values] = nonprompt_for_items()
                aborted = true;
                displayNames = cell(0,1);
                values = cell(0,1);
            end

            function [] = nop(varargin)
            end
        end
        
        
        function [btnSelectAll] = make_select_all_button_template()
            import Fancy.UI.FancyList.FancyListMgrBtn;
            btnSelectAll = FancyListMgrBtn('Select all', @(~, ~, flm) flm.select_all());
        end
        
        function [btnDeselectAll] = make_deselect_all_button_template()
            import Fancy.UI.FancyList.FancyListMgrBtn;
            btnDeselectAll = FancyListMgrBtn('Deselect all', @(~, ~, flm) flm.deselect_all());
        end
        
        function [btnAddItems] = make_add_items_button_template(fn_prompt_for_items, disableAddItemsButton)
            import Fancy.UI.FancyList.FancyListMgrBtn;
            if nargin < 2
                disableAddItemsButton = false;
            end
            
            function add_items_from_prompt(flm, fn_prompt_for_items)
                [aborted, itemDispTextCells, itemCells] = fn_prompt_for_items();
                if not(aborted)
                    flm.add_list_items(itemDispTextCells, itemCells);
                end
            end
            buttonText = 'Add item(s)';
            callback = @(~, ~, flm) add_items_from_prompt(flm, fn_prompt_for_items);
            if disableAddItemsButton
                btnAddItems = FancyListMgrBtn(buttonText, callback, 'off');
            else
                btnAddItems = FancyListMgrBtn(buttonText, callback);
            end
        end
        
        function [btnRemoveItems] = make_try_remove_items_button_template(warnBeforeRemovals, fn_on_items_removal)
            import Fancy.UI.FancyList.FancyListMgrBtn;

            function [] = try_remove_selected_items(flm, warnBeforeRemovals, fn_on_items_removal)
                indicesForRemoval = flm.get_selected_indices();


                itemsForRemoval = flm.get_display_value_pair_list(indicesForRemoval);

                if warnBeforeRemovals 
                    numDeletionItems = size(itemsForRemoval, 1);
                    if numDeletionItems == 0
                        msgbox('No entries were selected for removal!', 'Removal Failure');
                        return;
                    end
                    if numDeletionItems == 1
                        question = 'Are you sure you wish to remove the selected entry?';
                    else
                        question = sprintf('Are you sure you wish to remove the %d selected entries?', numDeletionItems);
                    end
                    optConfirmRemoval = 'Yes';
                    optAbortRemoval = 'No';
                    optDefault = optAbortRemoval;
                    choice = questdlg(question, 'Removal Confirmation', optConfirmRemoval, optAbortRemoval, optDefault);
                    continueRemoval = strcmp(choice, optConfirmRemoval);

                    if not(continueRemoval)
                        return;
                    end
                end
                fn_on_items_removal(itemsForRemoval);

                flm.remove_list_items(indicesForRemoval);
            end
            btnRemoveItems = FancyListMgrBtn('Remove item(s)', @(~, ~, flm) try_remove_selected_items(flm, warnBeforeRemovals, fn_on_items_removal));
        end
    end
end
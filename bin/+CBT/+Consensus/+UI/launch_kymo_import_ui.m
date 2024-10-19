function [lm] = launch_kymo_import_ui(hPanelKymoImport, tsCBC)
    % launch_kymo_import_ui -
    %   adds a tab with list management UI/functionality  for
    %    kymographs

    import Fancy.UI.FancyList.FancyListMgr;
    lm = FancyListMgr();
    lm.set_ui_parent(hPanelKymoImport);
    lm.make_ui_items_listbox();
    
    import Fancy.UI.FancyList.FancyListMgrBtnSet;
    flmbs1 = FancyListMgrBtnSet();
    flmbs1.NUM_BUTTON_COLS = 2;
    flmbs1.add_button(FancyListMgr.make_select_all_button_template());
    flmbs1.add_button(FancyListMgr.make_deselect_all_button_template());
    flmbs1.add_button(make_add_kymos_from_dbm_btn(tsCBC));
    flmbs1.add_button(make_add_kymos_directly_btn(tsCBC));

    flmbs2 = FancyListMgrBtnSet();
    flmbs2.NUM_BUTTON_COLS = 1;
    flmbs2.add_button(make_remove_kymos_btn());
    
    lm.add_button_sets(flmbs1, flmbs2);
    
    function [btnAddKymos] = make_add_kymos_directly_btn(ts)
        import Fancy.UI.FancyList.FancyListMgrBtn;
        btnAddKymos = FancyListMgrBtn(...
            'Add kymographs directly', ...
            @(~, ~, lm) on_add_kymos_directly(lm, ts));
        
        
        function [] = on_add_kymos_directly(lm, ts)

            import OldDBM.General.Import.import_raw_kymos;
            [rawKymos, rawKymoFilepaths] = import_raw_kymos();

            import OldDBM.General.Import.prompt_files_bps_per_pixel_wrapper;
            [pixelsWidths_bps] = prompt_files_bps_per_pixel_wrapper(rawKymoFilepaths, ts);

            % todo: generate kymo structs directly instead of using data
            %  wrapper
            
            import OldDBM.General.DataWrapper;
            dbmODW = DataWrapper();
            
            import OldDBM.General.Import.set_raw_kymo_data;
            set_raw_kymo_data(dbmODW, rawKymos, rawKymoFilepaths, pixelsWidths_bps);
            
            import OldDBM.General.Export.extract_kymo_structs;
            kymoStructs = extract_kymo_structs(dbmODW);
            
            
            unalignedKymos = cellfun(...
                @(kymoStruct) kymoStruct.unalignedKymo, ...
                kymoStructs, ...
                'UniformOutput', false);
            kymoNames = cellfun(...
                @(kymoStruct) kymoStruct.displayName, ...
                kymoStructs, ...
                'UniformOutput', false);
            
            hTabUnalignedKymos = get_unaligned_kymos_tab(ts);
            hPanelUnalignedKymos = uipanel('Parent', hTabUnalignedKymos);
            delete(allchild(hPanelUnalignedKymos));
            import OldDBM.Kymo.UI.show_kymos_in_grid;
            show_kymos_in_grid(hPanelUnalignedKymos, unalignedKymos, kymoNames);

            lm.add_list_items(kymoNames, kymoStructs);
        end
    end
    
    
    function [btnAddKymos] = make_add_kymos_from_dbm_btn(ts)
        import Fancy.UI.FancyList.FancyListMgrBtn;
        btnAddKymos = FancyListMgrBtn(...
            'Add kymographs from  DBM sessions', ...
            @(~, ~, lm) on_add_kymos_from_dbm(lm, ts));

        function [aborted, kymoNames, kymoStructs] = prompt_kymos_from_DBM_session()
            % prompt_kymos_from_DBM_session - get kymo structs
            %   extracted from formatted DBM session .mat files provided
            %   by a prompt to the user

            import CBT.Consensus.Import.prompt_dbm_session_filepath;
            [aborted, sessionFilepath] = prompt_dbm_session_filepath();

            if aborted
                kymoStructs = cell(0, 1);
                kymoNames = cell(0, 1);
                return;
            end

            import OldDBM.General.Import.try_loading_from_session_file;
            [dbmODW, dbmOSW] = try_loading_from_session_file(sessionFilepath);
            

            import OldDBM.General.Export.DataExporter;
            dbmDE = DataExporter(dbmODW, dbmOSW);
            kymoStructs = dbmDE.extract_kymo_structs();
            kymoNames = cellfun(...
                @(kymoStruct) kymoStruct.displayName, ...
                kymoStructs, ...
                'UniformOutput', false);
            
        end
        function [] = on_add_kymos_from_dbm(lm, ts)
            [aborted, kymoNames, kymoStructs] = prompt_kymos_from_DBM_session();
            if aborted
                return;
            end
            
            unalignedKymos = cellfun(...
                @(kymoStruct) kymoStruct.unalignedKymo, ...
                kymoStructs, ...
                'UniformOutput', false);
            
            hTabUnalignedKymos = get_unaligned_kymos_tab(ts);
            hPanelUnalignedKymos = uipanel('Parent', hTabUnalignedKymos);
            delete(allchild(hPanelUnalignedKymos));
            import OldDBM.Kymo.UI.show_kymos_in_grid;
            show_kymos_in_grid(hPanelUnalignedKymos, unalignedKymos, kymoNames);

            lm.add_list_items(kymoNames, kymoStructs);
        end
    end

    function [hTabUnalignedKymos] = get_unaligned_kymos_tab(tsCBC)
        persistent localhTabUnalignedKymos;
        if isempty(localhTabUnalignedKymos) || not(isvalid(localhTabUnalignedKymos))
            hTabUnalignedKymos = tsCBC.create_tab('Unaligned Kymos');
            localhTabUnalignedKymos = hTabUnalignedKymos;
        else
            hTabUnalignedKymos = localhTabUnalignedKymos;
        end
    end

    function [btnRemoveKymos] = make_remove_kymos_btn()
        import Fancy.UI.FancyList.FancyListMgrBtn;
        btnRemoveKymos = FancyListMgrBtn(...
            'Remove selected kymographs', ...
            @(~, ~, lm) on_remove_selected_kymos(lm));
        function [] = on_remove_selected_kymos(lm)
            lm.remove_selected_items();
        end
    end
end
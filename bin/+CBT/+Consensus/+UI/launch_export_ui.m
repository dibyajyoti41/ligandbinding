function [] = launch_export_ui(ts, exportDisplayNames, exportValues)
    % LAUNCH_EXPORT_UI
    %   adds a tab with list management UI/functionality for
    %   exportable data
    %  ts
    %   tabbed screen handle
    %  exportDisplayNames
    %   cell array of display names for values
    %  exportValues
    %   cell array of values

    tabTitle = 'Export Cluster Consensuses';
    hTab = ts.create_tab(tabTitle);
    hTabPanel = uipanel(hTab, 'Position', [0, 0, 1, 1]);
    ts.select_tab(hTab);

    import Fancy.UI.FancyList.FancyListMgrBtnSet;
    flmbs = FancyListMgrBtnSet();
    flmbs.NUM_BUTTON_COLS = 1;
    flmbs.add_button(make_export_cluster_mat_btn());
    flmbs.add_button(make_export_cluster_tiffs_btn());
    flmbs.add_button(make_export_cluster_tsvs_btn());
    flmbs.add_button(make_export_stretch_data_btn());

    function [btnExportClusterMat] = make_export_cluster_mat_btn()
        import Fancy.UI.FancyList.FancyListMgrBtn;
        btnExportClusterMat = FancyListMgrBtn(...
            'Export Selected Consensuses', ...
            @(~, ~, lm) on_export_selected_cluster_mat(lm));
        function [] = on_export_selected_cluster_mat(lm)
            import CBT.Consensus.Export.export_cluster_mat;

            selectedItems = lm.get_selected_list_items();
            selectedItemNames = selectedItems(:, 1);
            selectedItemValues = selectedItems(:, 2);
            numSelectedItems = length(selectedItemNames);

            for selectedItemNum=1:numSelectedItems
                clusterKey = selectedItemNames{selectedItemNum};
                clusterConsensusData = selectedItemValues{selectedItemNum};

                matFilepath = export_cluster_mat(clusterConsensusData, clusterKey);

                fprintf('Saved cluster consensus data for cluster ''%s'' to ''%s''\n', clusterKey, matFilepath);
            end
        end            
    end

    function [btnExportClusterTiffs] = make_export_cluster_tiffs_btn()
        import Fancy.UI.FancyList.FancyListMgrBtn;
        btnExportClusterTiffs = FancyListMgrBtn(...
            'Export Selected Consensuses'' Barcodes', ...
            @(~, ~, lm) on_export_selected_cluster_tiffs(lm));

        function [] = on_export_selected_cluster_tiffs(lm)
            import CBT.Consensus.Export.export_cluster_tiff;

            selectedItems = lm.get_selected_list_items();
            selectedItemNames = selectedItems(:, 1);
            selectedItemValues = selectedItems(:, 2);
            numSelectedItems = length(selectedItemNames);

            for selectedItemNum = 1:numSelectedItems
                clusterKey = selectedItemNames{selectedItemNum};
                clusterConsensusData = selectedItemValues{selectedItemNum};

                tiffFilepath = export_cluster_tiff(clusterConsensusData, clusterKey);

                fprintf('Saved cluster consensus image for cluster ''%s'' to ''%s''\n', clusterKey, tiffFilepath);
            end
        end

    end

    function [btnExportClusterTsvs] = make_export_cluster_tsvs_btn()
        import Fancy.UI.FancyList.FancyListMgrBtn;
        btnExportClusterTsvs = FancyListMgrBtn(...
            'Export Selected Consensuses'' TSV Spreadsheets', ...
            @(~, ~, lm) on_export_selected_cluster_tsvs(lm));

        function [] = on_export_selected_cluster_tsvs(lm)
            import CBT.Consensus.Export.export_selected_cluster_tsv;
            selectedItems = lm.get_selected_list_items();
            selectedItemNames = selectedItems(:, 1);
            selectedItemValues = selectedItems(:, 2);
            numSelectedItems = length(selectedItemNames);

            for selectedItemNum = 1:numSelectedItems
                clusterKey = selectedItemNames{selectedItemNum};
                clusterConsensusData = selectedItemValues{selectedItemNum};

                export_selected_cluster_tsv(clusterConsensusData, clusterKey);

                fprintf('Saved cluster consensus spreadsheet for cluster ''%s''\n', clusterKey);
            end
        end

    end

    function [btnExportStretchData] = make_export_stretch_data_btn()
        import Fancy.UI.FancyList.FancyListMgrBtn;
        btnExportStretchData = FancyListMgrBtn(...
            'Export Selected Consensuses'' Stretch Data TSV Spreadsheet', ...
            @(~, ~, lm) on_export_stretch_data_tsv(lm));

        function [] = on_export_stretch_data_tsv(lm)
            import CBT.Consensus.Export.export_stretch_data_tsv;

            selectedItems = lm.get_selected_list_items();
            selectedItemNames = selectedItems(:, 1);
            selectedItemValues = selectedItems(:, 2);
            numSelectedItems = length(selectedItemNames);

            for selectedItemNum = 1:numSelectedItems
                clusterKey = selectedItemNames{selectedItemNum};
                clusterConsensusData = selectedItemValues{selectedItemNum};

                export_stretch_data_tsv(clusterConsensusData, clusterKey);

            end
        end

    end

    import Fancy.UI.FancyList.FancyListMgr;
    lm = FancyListMgr();
    lm.set_ui_parent(hTabPanel);
    lm.make_ui_items_listbox();
    lm.add_list_items(exportDisplayNames, exportValues);
    lm.add_button_sets(flmbs);
end

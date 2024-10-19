function [] = add_consensus_btns_to_kymo_list_mgr(lm, ts, cache)
    if nargin < 3
        cache = containers.Map();
    end

    import Fancy.UI.FancyList.FancyListMgrBtnSet;
    flmbs2 = FancyListMgrBtnSet();
    flmbs2.NUM_BUTTON_COLS = 1;
    flmbs2.add_button(make_ensure_alignment_for_selected_kymos_btn(ts));
    flmbs2.add_button(make_generate_barcodes_for_selected_kymos_btn());
    flmbs2.add_button(make_consensus_btn(ts));

    lm.add_button_sets(flmbs2);

    function [hTabSelectedAlignedKymos] = get_aligned_selected_kymos_tab(ts)
        persistent localhTabSelectedAlignedKymos;
        if isempty(localhTabSelectedAlignedKymos) || not(isvalid(localhTabSelectedAlignedKymos))
            hTabSelectedAlignedKymos = ts.create_tab('Selected Aligned Kymos');
            localhTabSelectedAlignedKymos = hTabSelectedAlignedKymos;
        else
            hTabSelectedAlignedKymos = localhTabSelectedAlignedKymos;
        end
    end

    function [lm] = on_ensure_alignment_for_selected_kymos(lm, ts)
        import CBT.Consensus.Import.Helper.ensure_alignment_for_selected_kymos;
        [lm, kymoNames, alignedKymos] = ensure_alignment_for_selected_kymos(lm);

        hTabSelectedAlignedKymos = get_aligned_selected_kymos_tab(ts);
        hPanelAlignedKymos = uipanel('Parent', hTabSelectedAlignedKymos);
        delete(allchild(hPanelAlignedKymos));
        import OldDBM.Kymo.UI.show_kymos_in_grid;
        show_kymos_in_grid(hPanelAlignedKymos, alignedKymos, kymoNames);
    end

    function [btnEnsureAlignment] = make_ensure_alignment_for_selected_kymos_btn(ts)
        import Fancy.UI.FancyList.FancyListMgrBtn;
        btnEnsureAlignment = FancyListMgrBtn(...
            'Ensure Alignment of Selected Kymographs', ...
            @(~, ~, lm) on_ensure_alignment_for_selected_kymos(lm, ts));
    end

    function [btnEnsureAlignment] = make_generate_barcodes_for_selected_kymos_btn()
        import Fancy.UI.FancyList.FancyListMgrBtn;
        import CBT.Consensus.Import.Helper.generate_barcodes_for_selected_kymos;
        btnEnsureAlignment = FancyListMgrBtn(...
            'Generate Barcodes from Selected Kymographs', ...
            @(~, ~, lm) generate_barcodes_for_selected_kymos(lm));
    end

    function [btnGenerateConsensus] = make_consensus_btn(ts)
        function on_make_consensus(lm, ts)
            import CBT.Consensus.UI.Helper.generate_consensus_for_selected;
            [consensusStruct, cache] = generate_consensus_for_selected(lm, cache);
            
            import CBT.Consensus.Import.load_consensus_results;
            load_consensus_results(ts, consensusStruct);
        end

        import Fancy.UI.FancyList.FancyListMgrBtn;
        btnGenerateConsensus = FancyListMgrBtn(...
            'Generate Consensus from Selected Kymographs', ...
            @(~, ~, lm) on_make_consensus(lm, ts));
    end
end
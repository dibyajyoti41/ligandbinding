function [cache] = create_import_tab(hMenuParent, tsHCC, tabTitle,cache )
    % create_import_tab 
    if nargin < 4
        cache = containers.Map();
    end

    % create main tab for the analysis
    hTabKymoImport = tsHCC.create_tab(strcat([tabTitle ' import tab']));
    tsHCC.select_tab(hTabKymoImport);
    hPanelKymoImport = uipanel(hTabKymoImport);

    % import kymographs
    import Fancy.UI.Templates.launch_import_ui;
    [lm,cache] = launch_import_ui(hTabKymoImport,hPanelKymoImport,tsHCC,tabTitle,cache);
end


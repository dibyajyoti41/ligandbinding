function [ hMenuParent, tsHCA ] = create_figure_window( namefig, nametab )
    % creates figure window
    
    hFig = figure('Name', namefig, ...
        'Units', 'normalized', ...
        'OuterPosition', [0 0 1 1], ...
        'NumberTitle', 'off', ...
        'MenuBar', 'none' ...
    );

    hMenuParent = hFig;
    hPanel = uipanel('Parent', hFig);
    import Fancy.UI.FancyTabs.TabbedScreen;
    ts = TabbedScreen(hPanel);
    
    hTabHCA = ts.create_tab(nametab);
    ts.select_tab(hTabHCA);
    hPanelHCA = uipanel('Parent', hTabHCA);
    tsHCA = TabbedScreen(hPanelHCA);
    
    

end


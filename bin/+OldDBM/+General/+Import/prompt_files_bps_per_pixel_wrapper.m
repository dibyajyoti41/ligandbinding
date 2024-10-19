function [pixelsWidths_bps] = prompt_files_bps_per_pixel_wrapper(filepaths, ts)
    if nargin < 2
        ts = [];
    end
    if isempty(ts)
        hFigBpsPerPixel = figure('Name', 'Input Basepairs/Pixel');
        import Fancy.UI.FancyPositioning.try_maximize_figure;
        try_maximize_figure(hFigBpsPerPixel);
        hPanelBpsPerPixel = uipanel(hFigBpsPerPixel);
        import Fancy.UI.FancyTabs.TabbedScreen;
        ts = TabbedScreen(hPanelBpsPerPixel);
    end
    if isempty(filepaths)
        defaultVal = -1;
        pixelsWidths_bps = zeros(0, 1) + defaultVal;
    else
        hBpsPerPixelTab = ts.create_tab('bps/pixel');
        hBpsPerPixelPanel = uipanel(...
            'Parent', hBpsPerPixelTab, ...
            'Position', [0, 0, 1, 1]);
        ts.select_tab(hBpsPerPixelTab);
        import OldDBM.General.Import.prompt_files_bps_per_pixel;
        [pixelsWidths_bps, errorMsg] = prompt_files_bps_per_pixel(filepaths, hBpsPerPixelPanel);
        delete(hBpsPerPixelTab);
        waitfor(hBpsPerPixelTab);
        if not(isempty(errorMsg))
            disp(errorMsg);
            return;
        end
        pixelsWidths_bps = pixelsWidths_bps(:);
    end
end
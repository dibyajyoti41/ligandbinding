function [maximizationSuccessTF] = try_maximize_figure(hFig)
    % TRY_MAXIMIZE_FIGURE - try to maximize the figure
    %
    % Impressively, as of 2015, Matlab still provides no good,
    %   documented, or supported way of doing this, so we
    %   must resort to ugly, fragile, unreliable, or undocumented
    %   hacks with fallbacks
    % This is expected to stop working on some version. Hopefully Mathworks
    %   will provide some window maximize functionality by then
    %  https://mathworks.com/support/contact_us/dev/javaframe.html
    %
    % Inputs:
    %   hFig
    %    handle for a figure
    %
    % Outputs:
    %  maximizationSuccessTF
    %    true if the figure is successfully maximized
    %
    % Side-effects:
    %   if successful, maximizes the figure
    %
    % Authors:
    %   Saair Quaderi
    
    validateattributes(hFig, {'matlab.ui.Figure'}, {}, 1);
    
    maximizationSuccessTF = false;
    % Try getting undocumented Matlab internal feature (JavaFrame)
    %  which will *supposedly* be removed at some point]
    if isprop(hFig, 'JavaFrame')
        try
            drawnow;
            warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
            jFrame = get(hFig, 'JavaFrame');
            jFrame.setMaximized(true);
            maximizationSuccessTF = true;
        catch
            % No harm, no foul, right?
        end
    end
end

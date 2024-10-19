function [hAxes] = generate_axes_grid(hParent, numAxes)
    import Fancy.UI.FancyPositioning.FancyGrid.position_ui_elems_in_grid;
    if numAxes == 0
        hAxes = gobjects(0);
        return;
    end
    maxNumAxesPerRow = ceil(sqrt(numAxes));
    defaultArgStruct = struct('Parent', hParent, 'Position', [0, 0, 1, 1], 'Units', 'Normalized', 'XTick', [], 'YTick', [], 'XTickLabel', '', 'YTickLabel', '');
    overrideArgStructs = arrayfun(@(~) struct(), (1:numAxes)', 'UniformOutput', false);
    
    import Fancy.UI.generate_ui_elems;
    [hAxes] = generate_ui_elems(@axes, overrideArgStructs, defaultArgStruct);
    position_ui_elems_in_grid(hAxes, maxNumAxesPerRow);
end
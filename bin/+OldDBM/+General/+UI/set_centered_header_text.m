function [] = set_centered_header_text(hAxis, headerText, textColor, textBackgroundColor)
    if (nargin < 3) || isempty(textColor)
        textColor = [1 1 0];
    end
    if (nargin < 4) || isempty(textBackgroundColor)
        textBackgroundColor = 'none';
    end
    xLim = get(hAxis, 'XLim');
    yLim = get(hAxis, 'YLim');
    xCoord = ((xLim(2) - xLim(1))*.5);
    if strcmpi(get(hAxis, 'XDir'), 'reverse')
        xCoord = xLim(2) - xCoord;
    else
        xCoord = xLim(1) + xCoord;
    end
    yCoord =  ((yLim(2) - yLim(1))*.9);
    if strcmpi(get(hAxis, 'YDir'), 'reverse')
        yCoord = yLim(2) - yCoord;
    else
        yCoord = yLim(1) + yCoord;
    end

    axes(hAxis);
    text(  ...
        xCoord, yCoord, ...
        headerText, ...
        'Color', textColor, ...
        'BackgroundColor', textBackgroundColor, ...
        'HorizontalAlignment', 'center', ...
        'Interpreter', 'none', ...
        'FontWeight', 'bold', ...
        'FontSize', 11 ...
    );
end
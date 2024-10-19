function [figTitles, figNames, figNums, figHasNumTitleOnMask] = get_figure_titles(hFigs)
    % GET_FIGURE_TITLES - gets the titles of the figures
    %
    % Inputs:
    %   hFigs (optional):
    %     the handle of the figures (defaults to all available figures)
    %
    % Outputs:
    %   figTitles
    %     string cell array with the titles of the figure windows
    %   figNames
    %     string cell array with the names of the figures
    %   figNumbers
    %     string cell array with the numbers of the figures
    %   figHasNumTitleOn
    %     logical vector of whether each figure has NumberTitle set to 'on'
    %
    % Authors:
    %   Saair Quaderi
    if nargin < 1
        hFigs = findobj('Type', 'figure');
    else
        validateattributes(hFigs, {'matlab.ui.Figure'}, {}, 1);
    end
    
    figNames = arrayfun(@(f) get(f, 'Name'), hFigs, 'UniformOutput', false);
    figHasName = not(isempty(figNames));
    figNums = arrayfun(@(f) num2str(get(f, 'Number')), hFigs, 'UniformOutput', false);
    figHasNumTitleOnMask = arrayfun(@(f) strcmpi('on', get(f, 'NumberTitle')), hFigs);
    numFigs = numel(figNames);
    figTitles = cell(numFigs, 1);
    figNumTitles = cell(numFigs, 1);
    figNumTitles(figHasNumTitleOnMask) = strcat({'Figure '}, figNums(figHasNumTitleOnMask));
    figTitles(~figHasName & ~figHasNumTitleOnMask) = {''};
    figTitles(~figHasName & figHasNumTitleOnMask) = figNumTitles(~figHasName & figHasNumTitleOnMask);
    figTitles(figHasName & ~figHasNumTitleOnMask) = figNames(figHasName & ~figHasNumTitleOnMask);
    figTitles(figHasName & figHasNumTitleOnMask) = strcat(figNumTitles(figHasName & figHasNumTitleOnMask), {': '}, figNames(figHasName & figHasNumTitleOnMask));
end
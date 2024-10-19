function [] = close_tab_confirmation(hTab)
    % CLOSE_TAB_CONFIRMATION - deletes closest tab ancestor of
    %   the handle object
    %   presents the user the option to cancel
    %
    % Inputs:
    %   hTab
    %     a uitab handle
    %
    % Side-effects:
    %   closest tab ancestor of the handle object is deleted
    %    if user confirms this is intended through the popup
    %    menu
    %
    % Authors:
    %   Saair Quaderi
    
    validateattributes(hTab, {'matlab.ui.container.Tab'}, {'scalar'}, 1);
    tabTitle = get(hTab, 'Title');

    choice = menu('Are you sure you wish to close this tab?', 'Yes, close it', 'No, cancel');
    if choice == 1
        fprintf('Closing tab ''%s''\n', tabTitle);
        delete(hTab);
    end
end
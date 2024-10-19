function [hButtonContinue, hBpsPerPixelTable] = async_prompt_files_bps_per_pixel(filepaths, hPanel, onCompleteCallback)
    import Fancy.UI.FancyPositioning.get_pixel_height;

    function filename = get_filename(filepath)
        [~, filenameSansExt, fileExt] = fileparts(filepath);
        filename = [filenameSansExt, fileExt];
    end
    filenames = cellfun(@get_filename, filepaths, 'UniformOutput', false);
    numFiles = numel(filepaths);

    defaultFileBpsPerPixel = -1; % represents unknown value
    defaultFilesBpsPerPixel = repmat({defaultFileBpsPerPixel}, [numFiles, 1]);

    % Create the uitable
    defaultData = defaultFilesBpsPerPixel;
    rowNames = filenames;
    columnNames = {'<html>Basepairs/pixel (-1 = UNKNOWN)</html>'};
    columnFormats = {'numeric'};
    columnEditability = true;
    hBpsPerPixelTable = uitable(hPanel,...
        'RowName', rowNames,...
        'ColumnName', columnNames,...
        'ColumnFormat', columnFormats,...
        'ColumnEditable', columnEditability,...
        'Data', defaultData);
    hButtonContinue = uicontrol(hPanel, 'Style', 'pushbutton', 'String', 'Continue');

    PADDING_HEIGHT_PX = 20;
    buttonWidthPx = 200;
    buttonHeightPx = 30;
    % Set width and height
    tmpExtent = get(hBpsPerPixelTable, 'Extent');
    tableWidthPx = tmpExtent(3);
    tableHeightPx = tmpExtent(4);
    set(hBpsPerPixelTable, ...
        'Units', 'pixels', ...
        'Position', [0, buttonHeightPx + 2*PADDING_HEIGHT_PX, tableWidthPx, tableHeightPx]);

    set(hButtonContinue, ...
        'Units', 'pixels', ...
        'Position', [0, PADDING_HEIGHT_PX, buttonWidthPx, buttonHeightPx]);
    
%         %TODO: validate input changes
%         function on_cell_edit(~, ~)
%             % validate
%         end
%         iptaddcallback(paramTableHandle, 'CellEditCallback', {@on_cell_edit});


    incomplete = true;
    function continue_click(~, ~)
        if incomplete
            outData = get(hBpsPerPixelTable, 'Data');
            filesBpsPerPixel = cell2mat(outData(:,1));

            onCompleteCallback(filesBpsPerPixel);
            set(hBpsPerPixelTable, 'ColumnEditable', false);
            incomplete = false;
        end
    end
    iptaddcallback(hButtonContinue, 'Callback', {@continue_click});
end
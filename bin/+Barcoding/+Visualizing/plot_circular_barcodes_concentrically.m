function hGoSurface = plot_circular_barcodes_concentrically(hAxis, sanitizedBarcodes, settings)
    import Barcoding.Visualizing.colormap_kry;
    if (nargin < 3) || isempty(settings) || strcmpi(settings, 'default')
        settings = struct;
    end

    if not(iscell(sanitizedBarcodes))
        sanitizedBarcodes = {sanitizedBarcodes};
    end

    numBarcodes = length(sanitizedBarcodes);
    if numBarcodes < 1
        error('Must provide at least one barcode');
    end
    sanitizedBarcodes = cellfun(@(x) x(:)', sanitizedBarcodes, 'UniformOutput', false);



    if not(isfield(settings, 'BARCODE_WIDTH'))
        settings.BARCODE_WIDTH = 5;
    end
    if not(isfield(settings, 'GAP_WIDTH'))
        settings.GAP_WIDTH = settings.BARCODE_WIDTH/5;
    end
    if not(isfield(settings, 'FINAL_GAP_WIDTH'))
        settings.FINAL_GAP_WIDTH = 3*settings.GAP_WIDTH;
    end
    if not(isfield(settings, 'INNER_BAND_WIDTH'))
        settings.INNER_BAND_WIDTH = 8*settings.GAP_WIDTH;
    end
    if not(isfield(settings, 'BG_COLOR'))
%                 settings.BG_COLOR = [0.7, 0.75, 1.0]; 
        settings.BG_COLOR = [1.0, 1.0, 1.0];
    end

    if not(isfield(settings, 'COLORBAR_COLOR'))
        settings.COLORBAR_COLOR = [0 0 0 ];
    end
    if not(isfield(settings, 'COLOR_MAP'))
        settings.COLOR_MAP = colormap_kry();
    elseif ischar(settings.COLOR_MAP)
        try
            parenIdx = min(strfind(settings.COLOR_MAP,'('));
            if not(isempty(parenIdx))
                settings.COLOR_MAP = feval(settings.COLOR_MAP(1:parenIdx-1),str2double(settings.COLOR_MAP(parenIdx+1:end-1)));
            else
                settings.COLOR_MAP = feval(settings.COLOR_MAP);
            end
        catch
            settings.COLOR_MAP = [];
        end
    end

    widths = settings.GAP_WIDTH*ones(numBarcodes*2, 1);
    widths(1) = settings.INNER_BAND_WIDTH;
    widths(2:2:2*numBarcodes) = settings.BARCODE_WIDTH;
    widths(2*numBarcodes - 1) = settings.FINAL_GAP_WIDTH;
    radii = cumsum([0; widths]);


    barcodeLengths = cellfun(@(x) length(x), sanitizedBarcodes);
    barcodeLengths = barcodeLengths(:);
    maxLen = max(barcodeLengths);
    theta = -pi:(2*pi/maxLen):pi;
    Xs = radii*cos(theta);
    Ys = radii*sin(theta);
    Cs = NaN(2*numBarcodes + 1, maxLen + 1);

    for barcodeNum=1:numBarcodes
        currBarcode = [sanitizedBarcodes{barcodeNum}, NaN(1, maxLen - barcodeLengths(barcodeNum))];
        Cs(barcodeNum*2, :) = [currBarcode, currBarcode(1)];
    end
    hGoSurface = pcolor(hAxis, Xs, Ys, Cs);
    set(hAxis, 'color', settings.BG_COLOR, 'xtick', [], 'ytick', [], 'box', 'off');
    [~] = colorbar(hAxis, 'Color', settings.COLORBAR_COLOR, 'FontSize', 25);
    set(hAxis.Parent, 'BackgroundColor', settings.BG_COLOR);
    if not(isempty(settings.COLOR_MAP))
        set(get(hAxis, 'ColorSpace'), 'Colormap', settings.COLOR_MAP);
    end
    axis square;
    axis off;
    set(hAxis, 'XTickLabel', [], 'YTickLabel',[]); 
    grid off;
    shading flat; %shading interp;
end

function [edgeLength_px_unrounded, psfSigmaWidth_nm, deltaCut, pixelWidth_nm] = get_raw_px_edge_length(default_psfSigmaWidth_nm, default_deltaCut, default_pixelWidth_nm, skipOverridePrompt)
    if nargin < 4
        skipOverridePrompt = true;
    end
    psfSigmaWidth_nm = default_psfSigmaWidth_nm;
    deltaCut = default_deltaCut;
    pixelWidth_nm = default_pixelWidth_nm;
    if skipOverridePrompt
        edgeLength_px_unrounded = deltaCut * psfSigmaWidth_nm / pixelWidth_nm;
       return;
    end
    % Set a threshold for grouping barcodes.
    answer = inputdlg(...
        { ...
            'Point Spread Function Sigma Width (nm)', ...
            'Pixel width (nm)', ...
            'Delta Cut (edge length as multiple of PSF width in pixels))' ...
        }, ... % prompt
        'Bitmask Generation Parameters', ... % dialog title
        1, ... % number of lines
        { ...
            num2str(psfSigmaWidth_nm), ...
            num2str(pixelWidth_nm), ...
            num2str(deltaCut) ...
        }... % default value
    );
    if isempty(answer)
        edgeLength_px_unrounded = deltaCut * psfSigmaWidth_nm / pixelWidth_nm;
        return;
    else
        psfSigmaWidth_nm = str2double(answer{1});
        pixelWidth_nm = str2double(answer{2});
        deltaCut = str2double(answer{3});
    end
    isAcceptable = struct;
    isAcceptable.psfSigmaWidth_nm = (psfSigmaWidth_nm >= 0);
    isAcceptable.pixelWidth_nm = (pixelWidth_nm > 0);
    isAcceptable.deltaCut = (deltaCut >= 0);

    if all([isAcceptable.psfSigmaWidth_nm, isAcceptable.pixelWidth_nm, isAcceptable.deltaCut])
        edgeLength_px_unrounded = deltaCut * psfSigmaWidth_nm / pixelWidth_nm;
        return;
    end
    if not(isAcceptable.psfSigmaWidth_nm)
        psfSigmaWidth_nm = default_psfSigmaWidth_nm;
        warning('Bad input for PSF width! Try again!');
    end
    if not(isAcceptable.pixelWidth_nm)
        pixelWidth_nm = default_pixelWidth_nm;
        warning('Bad input pixel width! Try again!');
    end
    if not(isAcceptable.deltaCut)
        deltaCut = default_deltaCut;
        warning('Bad input for Delta Cut! Try again!');
    end
    edgeLength_px_unrounded = deltaCut * psfSigmaWidth_nm / pixelWidth_nm;
end
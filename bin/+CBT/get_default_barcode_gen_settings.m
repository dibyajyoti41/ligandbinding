function [defaultBarcodeGenSettings] = get_default_barcode_gen_settings()

    % Competitive binding parameters:
    defaultBarcodeGenSettings.concNetropsin_molar = 6; % Netropsin concentration, units molar
    defaultBarcodeGenSettings.concYOYO1_molar = 0.02; % YOYO-1 concentration, units molar
    defaultBarcodeGenSettings.concDNA = 2;
    % Micriscopy point spread function and pixel sampling related
    % parameters:
    defaultBarcodeGenSettings.meanBpExt_nm = 0.225; % mean extension length per basepair, units nm
    defaultBarcodeGenSettings.psfSigmaWidth_nm = 300; % width of std for gaussian approximation of point spread function, units nm
    defaultBarcodeGenSettings.pixelWidth_nm = 130; % width of pixels for the CCD camera, units nm
    
    % TODO: reconsider fundamental purpose of this parameter and whether
    %  it's variance meanBpExt_nm/uncertainty in psfSigmaWidth_nm/something
    %  else more specific that should be explicitly represented instead
    defaultBarcodeGenSettings.stretchFactor = 1; % factor for stretching barcode (1 = no stretching/compressing)
    
    % Relevant for experimental barcode bitmasks:
    defaultBarcodeGenSettings.deltaCut = 3; % multiple of point spread function standard deviation representing distance from experimental barcode cuts/edges where values are not considered sufficiently trustworthy

    
    defaultBarcodeGenSettings.isLinearTF = false; % whether to treat non-circularly (e.g. not using circular convolution of PSF)
    defaultBarcodeGenSettings.widthSigmasFromMean = 4; % number of psf widths to make the hSize for the gaussian kernel for PSF
    defaultBarcodeGenSettings.computeFreeConcentrations = 1; % whether to compute free concentrations
    defaultBarcodeGenSettings.model = 'literature'; % theory model to use
    % import CBT.Import.Helpers.read_CBT_settings_struct;
    % cbtSettingsStruct = read_CBT_settings_struct();
    % 
    % import Fancy.Utils.extract_fields;
    % [ ...
    %     defaultBarcodeGenSettings.concNetropsin_molar, ...
    %     defaultBarcodeGenSettings.concYOYO1_molar ...
    % ] = extract_fields(cbtSettingsStruct.cbtheory, { ...
    %     'NETROPSINconc', ...
    %     'YOYO1conc' ...
    %     });
    % 
    % [ ...
    %     defaultBarcodeGenSettings.psfSigmaWidth_nm, ...
    %     defaultBarcodeGenSettings.meanBpExt_nm, ...
    %     defaultBarcodeGenSettings.pixelWidth_nm ...
    % ] = extract_fields(cbtSettingsStruct.converttheory, { ...
    %     'PSF_width', ...
    %     'nmPerBps', ...
    %     'nmPerPixel' ...
    %     });
end
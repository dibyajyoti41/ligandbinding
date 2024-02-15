function [sets] = import_hca_settings(setsName)
   % import settings from file

    % settings
    sets = ini2struct( setsName);

    sets.bitmasking.untrustedPx = sets.bitmasking.deltaCut*sets.bitmasking.psfSigmaWidth_nm/sets.bitmasking.prestretchPixelWidth_nm;

    % default edge detection settings
    sets.edgeDetectionSettings = ini2struct( 'default_edge_detection.txt');

    import OptMap.MoleculeDetection.EdgeDetection.get_default_edge_detection_settings;
    sets.edgeDetectionSettings = get_default_edge_detection_settings(sets.skipDoubleTanhAdjustment);
    sets.edgeDetectionSettings.method = 'Otsu';
end


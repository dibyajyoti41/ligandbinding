function [experimentCurveStructs] = extract_experimental_structs(dbmODW)
    % extract_experimental_structs - get experiment curve structs
    %   built from data extracted from the DBM session mainStruct
    %   provided as an input
    import Fancy.Utils.data_hash;

    [fgKymoTimeAvgs, fgKymoTimeAvgFileIdxs, fgKymoTimeAvgFileMoleculeIdxs] = dbmODW.get_all_existing_fg_kymo_time_avgs();
    srcFilepaths = dbmODW.get_molecule_src_filepaths(fgKymoTimeAvgFileIdxs);
    srcFilenames = dbmODW.get_molecule_src_filenames(fgKymoTimeAvgFileIdxs);

    experimentDataType = 'Single Kymo Av';
    experimentCurves_pxRes = cellfun(@zscore, fgKymoTimeAvgs, 'UniformOutput', false);
    dataHashes = cellfun(@data_hash, experimentCurves_pxRes, 'UniformOutput', false);
    displayNames = strcat(srcFilenames, arrayfun(@(fileMoleculeIdx) [' Mol #', num2str(fileMoleculeIdx)], fgKymoTimeAvgFileMoleculeIdxs, 'UniformOutput', false));

    experimentCurveStructs = struct(...
        'filePath', srcFilepaths,...
        'displayName', displayNames,...
        'experimentCurve_pxRes', experimentCurves_pxRes,...
        'dataHash', dataHashes,...
        'type', experimentDataType...
        );
    experimentCurveStructs = arrayfun(@(experimentCurveStruct) experimentCurveStruct, experimentCurveStructs, 'UniformOutput', false);
end
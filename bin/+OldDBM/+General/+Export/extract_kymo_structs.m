function [kymoStructs] = extract_kymo_structs(dbmODW)
    % extract_kymo_structs - get kymo structs
    %   built from data extracted from the DBM session
    %   provided as an input
    
    [rawKymos, rawKymoFileIdxs, rawKymoFileMoleculeIdxs] = dbmODW.get_all_existing_raw_kymos();
    [alignedKymos, alignedKymosStretchFactors, shiftAlignedKymos] = dbmODW.get_aligned_kymos(rawKymoFileIdxs, rawKymoFileMoleculeIdxs);

    moleculeStatuses = dbmODW.get_molecule_statuses(rawKymoFileIdxs, rawKymoFileMoleculeIdxs);
    passesFilters = moleculeStatuses.passesFilters;
    if not(iscell(passesFilters))
        passesFilters = arrayfun(@(tf) tf, passesFilters, 'UniformOutput', false);
    end
    srcFilepaths = dbmODW.get_molecule_src_filepaths(rawKymoFileIdxs);
    srcFilenames = dbmODW.get_molecule_src_filenames(rawKymoFileIdxs);

    dataType = 'Kymo';
    import Fancy.Utils.data_hash;
    dataHashes = cellfun(@data_hash, rawKymos, 'UniformOutput', false);
    displayNames = strcat(srcFilenames, arrayfun(@(fileMoleculeIdx) [' Mol #', num2str(fileMoleculeIdx)], rawKymoFileMoleculeIdxs, 'UniformOutput', false));
    pixelsWidths_bp = dbmODW.get_molecule_src_pixel_widths_in_bps(rawKymoFileIdxs);
    if not(iscell(pixelsWidths_bp))
        pixelsWidths_bp = arrayfun(@(pixelsWidth_bp) pixelsWidth_bp, pixelsWidths_bp, 'UniformOutput', false);
    end

    kymoStructs = struct(...
        'passesFilters', passesFilters,...
        'filePath', srcFilepaths,...
        'displayName', displayNames,...
        'unalignedKymo', rawKymos,...
        'dataHash', dataHashes,...
        'alignedKymos', alignedKymos,...
        'alignedKymosStretchFactors', alignedKymosStretchFactors,...
        'shiftAlignedKymos', shiftAlignedKymos,...
        'bpsPerPx', pixelsWidths_bp,...
        'type', dataType...
    );
    kymoStructs = arrayfun(@(kymoStruct) kymoStruct, kymoStructs, 'UniformOutput', false);
end
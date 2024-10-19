function [] = export_raw_kymos(dbmODW, defaultOutputDirpath)
    [rawKymos, rawKymoFileIdxs, rawKymoFileMoleculeIdxs] = dbmODW.get_all_existing_raw_kymos();
    if isempty(rawKymos)
        fprintf('No raw kymographs were present to be exported\n');
        return;
    end

    outputDirpath = uigetdir(defaultOutputDirpath, 'Select Directory to Save Raw Kymo Files');

    if isequal(outputDirpath, 0)
        return;
    end

    [srcFilenames] = dbmODW.get_molecule_src_filenames(rawKymoFileIdxs);

    numRawKymos = length(rawKymos);
    outputKymoFilepaths = cell(numRawKymos, 1);
    for rawKymoNum = 1:numRawKymos
        [~, srcFilenameNoExt, ~] = fileparts(srcFilenames{rawKymoNum});
        fileMoleculeIdx = rawKymoFileMoleculeIdxs(rawKymoNum);
        outputKymoFilename = sprintf('%s_molecule_%d_kymograph.tif', srcFilenameNoExt, fileMoleculeIdx);
        outputKymoFilepath = fullfile(outputDirpath, outputKymoFilename);
        outputKymoFilepaths{rawKymoNum} = outputKymoFilepath;
    end
    
    cellfun(@(rawKymo, outputKymoFilepath)...
        imwrite(rawKymo, outputKymoFilepath, 'tif'),...
        rawKymos, outputKymoFilepaths);
end
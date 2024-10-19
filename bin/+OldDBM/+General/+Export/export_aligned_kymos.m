function [] = export_aligned_kymos(dbmODW, defaultOutputDirpath)
    [alignedKymos, ~, ~, alignedKymoFileIdxs, alignedKymoFileMoleculeIdxs] = dbmODW.get_all_existing_aligned_kymos();
    if isempty(alignedKymos)
        disp('No aligned kymographs were present to be exported')
        return;
    end

    outputDirpath = uigetdir(defaultOutputDirpath, 'Select Directory to Save Aligned Kymo Files');

    if isequal(outputDirpath, 0)
        return;
    end

    [srcFilenames] = dbmODW.get_molecule_src_filenames(alignedKymoFileIdxs);

    numAlignedKymos = length(alignedKymos);
    for alignedKymoNum = 1:numAlignedKymos
        [~, srcFilename, ~] = fileparts(srcFilenames{alignedKymoNum});
        alignedKymo = alignedKymos{alignedKymoNum};
        fileMoleculeIdx = alignedKymoFileMoleculeIdxs(alignedKymoNum);
        outputAlignedKymoFilepath = fullfile(outputDirpath, [srcFilename '_molecule_' num2str(fileMoleculeIdx) '_aligned_kymograph.tif']);
        imwrite(alignedKymo, outputAlignedKymoFilepath, 'tif');
    end
end
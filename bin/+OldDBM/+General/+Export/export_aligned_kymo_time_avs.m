function [] = export_aligned_kymo_time_avs(dbmODW, defaultOutputDirpath)
    [alignedKymos, ~, ~, alignedKymoFileIdxs, alignedKymoFileMoleculeIdxs] = dbmODW.get_all_existing_aligned_kymos();
    if isempty(alignedKymos)
        fprintf('No aligned kymographs were present to have their time-averages exported\n')
        return;
    end

    outputDirpath = uigetdir(defaultOutputDirpath, 'Select Directory to Save Kymo Time-Average Text Files');

    if isequal(outputDirpath, 0)
        return;
    end

    [srcFilenames] = dbmODW.get_molecule_src_filenames(alignedKymoFileIdxs);

    numAlignedKymos = length(alignedKymos);
    timeAvgOutputFilepaths = cell(numAlignedKymos, 1);
    for alignedKymoNum = 1:numAlignedKymos
        [~, srcFilename, ~] = fileparts(srcFilenames{alignedKymoNum});
        fileMoleculeIdx = alignedKymoFileMoleculeIdxs(alignedKymoNum);
        timeAvgOutputFilename = sprintf('%s_molecule_%d_aligned_kymo_time_avg.txt', srcFilename, num2str(fileMoleculeIdx));
        timeAvgOutputFilepath = fullfile(outputDirpath, timeAvgOutputFilename);
        timeAvgOutputFilepaths{alignedKymoNum} = timeAvgOutputFilepath;
    end
    for alignedKymoNum = 1:numAlignedKymos
        alignedKymo = alignedKymos{alignedKymoNum};
        alignedKymoTimeAvg = mean(alignedKymo, 1);
        numPixels = length(alignedKymoTimeAvg);
        timeAvgOutputFilepath = timeAvgOutputFilepaths{alignedKymoNum};
        fid = fopen(timeAvgOutputFilepath, 'wt');
        for pixelNum=1:numPixels
            fprintf(fid, '%s ', num2str(alignedKymoTimeAvg(pixelNum)));
        end
        fclose(fid);
    end
end
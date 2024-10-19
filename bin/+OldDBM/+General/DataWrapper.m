classdef DataWrapper < handle
    %DATAWRAPPER Functions for maintaining compatibility with old
    %structure of DBM data
    
    properties %(Access = private)
        DBMMainstruct
        SessionFilepath
    end  
    methods
        function [dbmODW] = DataWrapper(mainStruct, sessionFilepath)
            if nargin < 1
                mainStruct = struct();
                mainStruct.fileMoleculeCell = {};
                mainStruct.fileCell = {};
            end
            if nargin < 2
                sessionFilepath = '';
            end
            
            if not(isfield(mainStruct, 'fileMoleculeCell')) || not(isfield(mainStruct, 'fileCell'))
                error('The main struct for the session is not in a recognized format');
            end
            
            if isfield(mainStruct, 'axesArr')
                mainStruct = rmfield(mainStruct, 'axesArr');
            end
            if isfield(mainStruct, 'fig')
                mainStruct = rmfield(mainStruct, 'fig');
            end
            
            dbmODW.DBMMainstruct = mainStruct;
            dbmODW.SessionFilepath = sessionFilepath;
        end
        
        
        function [] = update_data(dbmODW, dbmODW2)
            dbmODW.DBMMainstruct = dbmODW2.DBMMainstruct;
            dbmODW.SessionFilepath = dbmODW2.SessionFilepath;
        end
        
        function [] = set_filter_settings(dbmODW, filterSettings)
            if isfield(filterSettings, 'minInfoScoreThresh')
                dbmODW.DBMMainstruct.filterInfo = filterSettings.minInfoScoreThresh;
            end
            
            if isfield(filterSettings, 'minMoleculeLengthEstimateThresh_px')
                dbmODW.DBMMainstruct.filterLength = filterSettings.minMoleculeLengthEstimateThresh_px;
            end
        end


        function [filterSettings] = get_filter_settings(dbmODW)
            filterSettings = struct();
            
            minInfoScoreThresh = 0;
            if isfield(dbmODW.DBMMainstruct, 'filterInfo')
                minInfoScoreThresh = dbmODW.DBMMainstruct.filterInfo;
            end
            
            minMoleculeLengthEstimateThresh_px = 0;
            if isfield(dbmODW.DBMMainstruct, 'minMoleculeLengthEstimateThresh_px')
               minMoleculeLengthEstimateThresh_px = dbmODW.DBMMainstruct.filterLength;
            end
            
            filterSettings.minInfoScoreThresh = minInfoScoreThresh;
            filterSettings.minMoleculeLengthEstimateThresh_px = minMoleculeLengthEstimateThresh_px;
        end
        
        
        
        

        function [averagingWindowWidth] = get_averaging_window_width(dbmODW)
            averagingWindowWidth = 3; % TODO: avoid hardcoding this value
            if isfield(dbmODW.DBMMainstruct, 'windowWidth')
                averagingWindowWidth = dbmODW.DBMMainstruct.windowWidth;
            end
        end

        
        
        
        
        
        
        
        function [fileIdxs, fileMoleculeIdxs] = get_molecule_idxs(dbmODW)
            numFilesMolecules = cellfun(@(fileMolecules) numel(fileMolecules), dbmODW.DBMMainstruct.fileMoleculeCell);
            numMolecules = sum(numFilesMolecules);
            fileIdxs = NaN(numMolecules, 1);
            fileMoleculeIdxs = NaN(numMolecules, 1);
            moleculeNum = 0;
            numFiles = numel(numFilesMolecules);
            for fileNum = 1:numFiles
                numFileMolecules = numFilesMolecules(fileNum);
                currFileMoleculeIdxs = (1:numFileMolecules)';
                idxIdxs = moleculeNum + currFileMoleculeIdxs;
                fileIdxs(idxIdxs) = fileNum;
                fileMoleculeIdxs(idxIdxs) = currFileMoleculeIdxs;
                moleculeNum = moleculeNum + numFileMolecules;
            end
        end
        
        function [moleculeStatuses] = get_molecule_statuses(dbmODW, fileIdxs, fileMoleculeIdxs)
            [moleculeStructs] = dbmODW.get_molecule_structs(fileIdxs, fileMoleculeIdxs);
            numMolecules = numel(fileIdxs);
            moleculeStatuses = struct();
            moleculeStatuses.passesFilters = false(numMolecules, 1);
            moleculeStatuses.hasMovie = false(numMolecules, 1);
            moleculeStatuses.hasRawKymo = false(numMolecules, 1);
            moleculeStatuses.hasAlignedKymo = false(numMolecules, 1);
            moleculeStatuses.hasFgKymoTimeAvg = false(numMolecules, 1);
            for moleculeNum=1:numMolecules
                moleculeStruct = moleculeStructs{moleculeNum};

                if moleculeStruct.passesFilters
                    moleculeStatuses.passesFilters(moleculeNum) = true;
                end
                if isfield(moleculeStruct, 'frames') && not(isempty(moleculeStruct.frames))
                    moleculeStatuses.hasMovie(moleculeNum) = true;
                end
                if isfield(moleculeStruct, 'kymograph') && not(isempty(moleculeStruct.kymograph))
                    moleculeStatuses.hasRawKymo(moleculeNum) = true;
                end
                if isfield(moleculeStruct, 'alignedKymograph') && not(isempty(moleculeStruct.alignedKymograph))
                    moleculeStatuses.hasAlignedKymo(moleculeNum) = true;
                end
                if isfield(moleculeStruct, 'timeAvKymo') && not(isempty(moleculeStruct.timeAvKymo))
                    moleculeStatuses.hasFgKymoTimeAvg(moleculeNum) = true;
                end
            end
        end
        
        function [srcFilepath] = get_molecule_src_filepath(dbmODW, fileIdx)
            srcFilepath = dbmODW.DBMMainstruct.fileCell{fileIdx}.fileName;
        end
        
        function [srcFilename] = get_molecule_src_filename(dbmODW, fileIdx)
            [srcFilepath] = dbmODW.get_molecule_src_filepath(fileIdx);
            [~, name, ext] = fileparts(srcFilepath);
            srcFilename = [name, ext];
        end

        function [moleculeRectPosition] = get_molecule_rect_position(dbmODW, fileIdx, fileMoleculeIdx)
            moleculeRectPosition = [NaN, NaN, NaN, NaN];
            try
                windowWidth = dbmODW.get_averaging_window_width();
                colStartIdx = dbmODW.DBMMainstruct.fileCell{fileIdx}.regions(fileMoleculeIdx, 1);
                colEndIdx = dbmODW.DBMMainstruct.fileCell{fileIdx}.regions(fileMoleculeIdx, 2);
                rowCenterIdx = dbmODW.DBMMainstruct.fileCell{fileIdx}.locs(fileMoleculeIdx);
                rowStartIdx = rowCenterIdx - floor(windowWidth/2);
                rowEndIdx = rowCenterIdx + floor(windowWidth/2);
                moleculeRectPosition = [colStartIdx, rowStartIdx, colEndIdx - colStartIdx, rowEndIdx - rowStartIdx];
            catch
            end
        end
        
        function [moleculeRectPositions] = get_molecule_rect_positions(dbmODW, fileIdxs, fileMoleculeIdxs)
            moleculeRectPositions = arrayfun(@(fileIdx, fileMoleculeIdx) dbmODW.get_molecule_rect_position(fileIdx, fileMoleculeIdx), fileIdxs, fileMoleculeIdxs, 'UniformOutput', false);
        end
        
        function [srcFilepaths] = get_molecule_src_filepaths(dbmODW, fileIdxs)
            srcFilepaths = arrayfun(@(fileIdx) dbmODW.get_molecule_src_filepath(fileIdx), fileIdxs, 'UniformOutput', false);
        end

        function [srcFilenames] = get_molecule_src_filenames(dbmODW, fileIdxs)
            srcFilenames = arrayfun(@(fileIdx) dbmODW.get_molecule_src_filename(fileIdx), fileIdxs, 'UniformOutput', false);
        end
        
        
        function [pixelsWidth_bp] = get_molecule_src_pixel_width_in_bps(dbmODW, fileIdx)
            fileStruct = dbmODW.DBMMainstruct.fileCell{fileIdx};
            pixelsWidth_bp = NaN;
            if isfield(fileStruct, 'expBpsPerPixel')
                pixelsWidth_bp = fileStruct.expBpsPerPixel;
            end
        end
        
        function [] = set_molecule_src_pixel_width_in_bps(dbmODW, fileIdx, pixelsWidth_bp)
            fileStruct = dbmODW.DBMMainstruct.fileCell{fileIdx};
            fileStruct.expBpsPerPixel = pixelsWidth_bp;
        end
        
        
        function [pixelsWidths_bp] = get_molecule_src_pixel_widths_in_bps(dbmODW, fileIdxs)
            pixelsWidths_bp = arrayfun(@(fileIdx) dbmODW.get_molecule_src_pixel_width_in_bps(fileIdx), fileIdxs);
        end
        
        function [] = set_molecule_src_pixel_widths_in_bps(dbmODW, fileIdxs, pixelsWidths_bps)
            arrayfun(@(fileIdx, pixelsWidth_bp) dbmODW.set_molecule_src_pixel_width_in_bps(fileIdx, pixelsWidth_bp), fileIdxs, pixelsWidths_bps);
        end
        
        
        function [movieFrameTimeAvg] = get_molecule_src_movie_frame_time_avg(dbmODW, fileIdx)
            fileStruct = dbmODW.DBMMainstruct.fileCell{fileIdx};
            movieFrameTimeAvg = NaN(0, 0);
            if isfield(fileStruct, 'averagedImg')
                movieFrameTimeAvg = (fileStruct.averagedImg)';
            end
        end
        
        
        function [moleculeStruct] = get_molecule_struct(dbmODW, fileIdx, fileMoleculeIdx)
            moleculeStruct = dbmODW.DBMMainstruct.fileMoleculeCell{fileIdx}{fileMoleculeIdx};
        end
        
        function [] = set_molecule_struct(dbmODW, fileIdx, fileMoleculeIdx, moleculeStruct)
            dbmODW.DBMMainstruct.fileMoleculeCell{fileIdx}{fileMoleculeIdx} = moleculeStruct;
        end
        
        function [moleculeStructs] = get_molecule_structs(dbmODW, fileIdxs, fileMoleculeIdxs)
            numMolecules = numel(fileIdxs);
            moleculeStructs = cell(numMolecules, 1);
            for moleculeNum = 1:numMolecules
                fileIdx = fileIdxs(moleculeNum);
                fileMoleculeIdx = fileMoleculeIdxs(moleculeNum);
                moleculeStruct = dbmODW.get_molecule_struct(fileIdx, fileMoleculeIdx);
                moleculeStructs{moleculeNum} = moleculeStruct;
            end
        end
        
        function [] = set_molecule_structs(dbmODW, fileIdxs, fileMoleculeIdxs, moleculeStructs)
            numMolecules = numel(fileIdxs);
            for moleculeNum = 1:numMolecules
                fileIdx = fileIdxs(moleculeNum);
                fileMoleculeIdx = fileMoleculeIdxs(moleculeNum);
                moleculeStruct = moleculeStructs{moleculeNum};
                
                dbmODW.set_molecule_struct(fileIdx, fileMoleculeIdx, moleculeStruct);
            end
        end
        
        
        
        
        function [] = verify_thresholds(dbmODW)
            % VERIFY_THRESHOLDS - Goes through all the kymographs and picks 
            % ones that pass the thresholds set in setoptions
            %
            % Authors:
            %   Charleston Noble
            %   Saair Quaderi


            filterSettings = dbmODW.get_filter_settings();
            minInfoScoreThresh = filterSettings.minInfoScoreThresh;
            minMoleculeLengthEstimateThresh_px = filterSettings.minMoleculeLengthEstimateThresh_px;

            [fileIdxs, fileMoleculeIdxs] = dbmODW.get_molecule_idxs();
            numMolecules = length(fileMoleculeIdxs);
            infoScores = dbmODW.get_info_scores(fileIdxs, fileMoleculeIdxs);
            passesInfoScoreFilterVect = isnan(infoScores) | (infoScores > minInfoScoreThresh);


            % SEEMS THAT THIS ALWAYS RETURNS NAN's ? TODO: make sure that length
            % estimate is always computed

            passesMoleculeLengthFilterVect = true(numMolecules, 1);
            for moleculeNum = 1:numMolecules
                fileIdx = fileIdxs(moleculeNum);
                fileMoleculeIdx = fileMoleculeIdxs(moleculeNum);

                moleculeLengthEstimate_px = dbmODW.get_molecule_length(fileIdx, fileMoleculeIdx);

                passesMoleculeLengthFilter = isnan(moleculeLengthEstimate_px) | (moleculeLengthEstimate_px < minMoleculeLengthEstimateThresh_px);

                passesMoleculeLengthFilterVect(moleculeNum) = passesMoleculeLengthFilter;
            end

            passesFiltersVect = passesInfoScoreFilterVect & passesMoleculeLengthFilterVect;


            [moleculeStructs] = dbmODW.get_molecule_structs(fileIdxs, fileMoleculeIdxs);

            for moleculeNum = 1:numMolecules
                moleculeStruct = moleculeStructs{moleculeNum};
                moleculeStruct.passesFilters = passesFiltersVect(moleculeNum);
                moleculeStructs{moleculeNum} = moleculeStruct;
            end
            dbmODW.set_molecule_structs(fileIdxs, fileMoleculeIdxs, moleculeStructs);
        end
        
        function [] = update_filter_settings(dbmODW)
            filterSettings = dbmODW.get_filter_settings();
            defaultMinMoleculeLengthEstimateThresh_px = filterSettings.minMoleculeLengthEstimateThresh_px;
            defaultMinInfoScoreThresh = filterSettings.minInfoScoreThresh;

            minMoleculeLengthEstimateThresh_px = defaultMinMoleculeLengthEstimateThresh_px;
            minInfoScoreThresh = defaultMinInfoScoreThresh;

            molecularLengthThresholdPrompt = 'Molecule length threshold';
            moleculeKymoInfoScoreThresholdPrompt = 'Molecule information threshold';

            prompts = {molecularLengthThresholdPrompt; moleculeKymoInfoScoreThresholdPrompt};
            dialogTitle = 'Set filter options';
            numLines = 1;
            defaultValStrs = {num2str(defaultMinMoleculeLengthEstimateThresh_px); num2str(minInfoScoreThresh)};
            answers = inputdlg(prompts, dialogTitle, numLines, defaultValStrs);

            if isempty(answers)
                return;
            end

            molecularLengthThresholdAns = str2double(answers{1});
            moleculeKymoInfoScoreThresholdAns = str2double(answers{2});

            if isnan(molecularLengthThresholdAns)
                disp('Invalid molecular length threshold');
            else
                minMoleculeLengthEstimateThresh_px = molecularLengthThresholdAns;
            end
            if isnan(moleculeKymoInfoScoreThresholdAns)
                disp('Invalid info score threshold');
            else
                minInfoScoreThresh = moleculeKymoInfoScoreThresholdAns;
            end


            filterSettings = struct();
            filterSettings.minMoleculeLengthEstimateThresh_px = minMoleculeLengthEstimateThresh_px;
            filterSettings.minInfoScoreThresh = minInfoScoreThresh;
            dbmODW.set_filter_settings(filterSettings);
        end
        
        function [] = update_molecule_stats(dbmODW, fileIdx, fileMoleculeIdx, kymoStatsStruct, useMedianNotMeanForLengthStats)
            moleculeStruct = dbmODW.get_molecule_struct(fileIdx, fileMoleculeIdx);
            if useMedianNotMeanForLengthStats
                lengthEst = kymoStatsStruct.medianOfFramewiseMoleculeExts;
                lengthDevEst = kymoStatsStruct.madEstStdOfFramewiseMoleculeExts;
            else
                lengthEst = kymoStatsStruct.meanOfFramewiseMoleculeExts;
                lengthDevEst = kymoStatsStruct.stdOfFramewiseMoleculeExts;
            end
            moleculeStructOverrides.length = lengthEst;
            moleculeStructOverrides.lengthSTD = lengthDevEst;
            moleculeStructOverrides.leftEndArr = kymoStatsStruct.moleculeLeftEdgeIdxs;
            moleculeStructOverrides.rightEndArr = kymoStatsStruct.moleculeRightEdgeIdxs;
            moleculeStructOverrides.COM = kymoStatsStruct.meanUnroundedCenterOfMassIdx;
            moleculeStructOverrides.INT = kymoStatsStruct.meanFramewiseMoleculeIntensity;
            moleculeStructOverrides.stdINT = kymoStatsStruct.stdFramewiseMoleculeIntensity ;
            moleculeStructOverrides.BGint = kymoStatsStruct.meanNonMainMoleculePixelIntensity; % todo: fix if possible -- warning: could include intensities of other non-main molecule foreground in addition to background!

            import Fancy.Utils.merge_structs;
            moleculeStruct = merge_structs(moleculeStruct, moleculeStructOverrides);
            dbmODW.set_molecule_struct(fileIdx, fileMoleculeIdx, moleculeStruct);
        end
        
        function [] = update_molecules_stats(dbmODW, kymoStatsTable, useMedianNotMeanForLengthStats)
            rowfun( ...
                @(fileIdx, fileMoleculeIdx, kymoStatsStruct) ...
                    dbmODW.update_molecule_stats(fileIdx, fileMoleculeIdx, kymoStatsStruct, useMedianNotMeanForLengthStats), ...
                 kymoStatsTable);
        end
        
        
        function [miniMovie] = get_mini_movie(dbmODW, fileIdx, fileMoleculeIdx)
            miniMovie = [];
            [moleculeStruct] = dbmODW.get_molecule_struct(fileIdx, fileMoleculeIdx);
            if isfield(moleculeStruct, 'frames')
                miniMovie = moleculeStruct.frames;
            end
        end
        
        function [] = set_movie(dbmODW, fileIdx, fileMoleculeIdx, miniMovie)
            dbmODW.DBMMainstruct.fileMoleculeCell{fileIdx}{fileMoleculeIdx}.frames = miniMovie;
        end

        function [miniMovies] = get_mini_movies(dbmODW, fileIdxs, fileMoleculeIdxs)
            numMolecules = numel(fileIdxs);
            miniMovies = cell(numMolecules, 1);
            for moleculeNum=1:numMolecules
                fileIdx = fileIdxs(moleculeNum);
                fileMoleculeIdx = fileMoleculeIdxs(moleculeNum);
                miniMovie = dbmODW.get_mini_movie(fileIdx, fileMoleculeIdx);
                miniMovies{moleculeNum} = miniMovie;
            end
        end

        function [miniMovies, movieFileIdxs, movieFileMoleculeIdxs] = get_all_existing_mini_movies(dbmODW)
            [fileIdxs, fileMoleculeIdxs] = dbmODW.get_molecule_idxs();
            [moleculeStatuses] = dbmODW.get_molecule_statuses(fileIdxs, fileMoleculeIdxs);
            movieFileIdxs = fileIdxs(moleculeStatuses.hasMovie);
            movieFileMoleculeIdxs = fileMoleculeIdxs(moleculeStatuses.hasRawKymo);
            [miniMovies] = dbmODW.get_mini_movies(movieFileIdxs, movieFileMoleculeIdxs);
        end

        
        
        function [rawKymo] = get_raw_kymo(dbmODW, fileIdx, fileMoleculeIdx)
            rawKymo = [];
            [moleculeStruct] = dbmODW.get_molecule_struct(fileIdx, fileMoleculeIdx);
            if isfield(moleculeStruct, 'kymograph')
                rawKymo = moleculeStruct.kymograph;
            end
        end
        
        function [] = set_raw_kymo(dbmODW, fileIdx, fileMoleculeIdx, rawKymo)
            moleculeStruct = dbmODW.get_molecule_struct(fileIdx, fileMoleculeIdx);
            moleculeStruct.kymograph = rawKymo;
            dbmODW.set_molecule_struct(fileIdx, fileMoleculeIdx, moleculeStruct);
        end

        function [rawKymos] = get_raw_kymos(dbmODW, fileIdxs, fileMoleculeIdxs)
            numMolecules = numel(fileIdxs);
            rawKymos = cell(numMolecules, 1);
            for moleculeNum=1:numMolecules
                fileIdx = fileIdxs(moleculeNum);
                fileMoleculeIdx = fileMoleculeIdxs(moleculeNum);
                rawKymo = dbmODW.get_raw_kymo(fileIdx, fileMoleculeIdx);
                rawKymos{moleculeNum} = rawKymo;
            end
        end
        
        function [rawKymos, rawKymoFileIdxs, rawKymoFileMoleculeIdxs] = get_all_existing_raw_kymos(dbmODW)
            [fileIdxs, fileMoleculeIdxs] = dbmODW.get_molecule_idxs();
            [moleculeStatuses] = dbmODW.get_molecule_statuses(fileIdxs, fileMoleculeIdxs);
            rawKymoFileIdxs = fileIdxs(moleculeStatuses.hasRawKymo);
            rawKymoFileMoleculeIdxs = fileMoleculeIdxs(moleculeStatuses.hasRawKymo);
            [rawKymos] = dbmODW.get_raw_kymos(rawKymoFileIdxs, rawKymoFileMoleculeIdxs);
        end

        function [] = set_raw_kymos(dbmODW, fileIdxs, fileMoleculeIdxs, rawKymos)
            numKymos = numel(fileIdxs);
            for kymoNum=1:numKymos
                fileIdx = fileIdxs(kymoNum);
                fileMoleculeIdx = fileMoleculeIdxs(kymoNum);
                rawKymo = rawKymos{kymoNum};
                dbmODW.set_raw_kymo(fileIdx, fileMoleculeIdx, rawKymo);
            end
        end

        
       
        function [alignedKymo, alignedKymoStretchFactors, shiftAlignedKymo] = get_aligned_kymo(dbmODW, fileIdx, fileMoleculeIdx)
            alignedKymo = [];
            alignedKymoStretchFactors = [];
            shiftAlignedKymo = [];
            [moleculeStruct] = dbmODW.get_molecule_struct(fileIdx, fileMoleculeIdx);
            if isfield(moleculeStruct, 'alignedKymograph') && not(isempty(moleculeStruct.alignedKymograph))
                alignedKymo = moleculeStruct.alignedKymograph;
            end
            if isfield(moleculeStruct, 'alignedKymographStretchFactors') && not(isempty(moleculeStruct.alignedKymographStretchFactors))
                alignedKymoStretchFactors = moleculeStruct.alignedKymographStretchFactors;
            end
            if isfield(moleculeStruct, 'shiftAlignedKymograph') && not(isempty(moleculeStruct.shiftAlignedKymograph))
                shiftAlignedKymo = moleculeStruct.shiftAlignedKymograph;
            end
        end

        function [] = set_aligned_kymo(dbmODW, fileIdx, fileMoleculeIdx, alignedKymo, stretchFactorsMat, shiftAlignedKymo)
            moleculeStruct = dbmODW.get_molecule_struct(fileIdx, fileMoleculeIdx);
            moleculeStruct.alignedKymograph = alignedKymo;
            moleculeStruct.alignedKymographStretchFactors = stretchFactorsMat;
            moleculeStruct.shiftAlignedKymograph = shiftAlignedKymo;
            dbmODW.set_molecule_struct(fileIdx, fileMoleculeIdx, moleculeStruct); 
        end
        
        function [alignedKymos, alignedKymosStretchFactors, shiftAlignedKymos] = get_aligned_kymos(dbmODW, fileIdxs, fileMoleculeIdxs)
            numMolecules = numel(fileIdxs);
            alignedKymos = cell(numMolecules, 1);
            alignedKymosStretchFactors = cell(numMolecules, 1);
            shiftAlignedKymos = cell(numMolecules, 1);
            for moleculeNum = 1:numMolecules
                fileIdx = fileIdxs(moleculeNum);
                fileMoleculeIdx = fileMoleculeIdxs(moleculeNum);
                [alignedKymo, alignedKymoStretchFactors, shiftAlignedKymo] = dbmODW.get_aligned_kymo(fileIdx, fileMoleculeIdx);
                alignedKymos{moleculeNum} = alignedKymo;
                alignedKymosStretchFactors{moleculeNum} = alignedKymoStretchFactors;
                shiftAlignedKymos{moleculeNum} = shiftAlignedKymo;
            end
        end

        function [] = set_aligned_kymos(dbmODW, fileIdxs, fileMoleculeIdxs, alignedKymos, stretchFactorsMats, shiftAlignedKymos)
            numKymos = numel(fileIdxs);
            for kymoNum = 1:numKymos
                fileIdx = fileIdxs(kymoNum);
                fileMoleculeIdx = fileMoleculeIdxs(kymoNum);
                alignedKymo = alignedKymos{kymoNum};
                stretchFactorsMat = stretchFactorsMats{kymoNum};
                shiftAlignedKymo = shiftAlignedKymos{kymoNum};
                dbmODW.set_aligned_kymo(fileIdx, fileMoleculeIdx, alignedKymo, stretchFactorsMat, shiftAlignedKymo);
            end
        end

        function [alignedKymos, alignedKymosStretchFactors, shiftAlignedKymos, alignedKymoFileIdxs, alignedKymoFileMoleculeIdxs] = get_all_existing_aligned_kymos(dbmODW)
            [fileIdxs, fileMoleculeIdxs] = dbmODW.get_molecule_idxs();
            [moleculeStatuses] = dbmODW.get_molecule_statuses(fileIdxs, fileMoleculeIdxs);
            alignedKymoFileIdxs = fileIdxs(moleculeStatuses.hasAlignedKymo);
            alignedKymoFileMoleculeIdxs = fileMoleculeIdxs(moleculeStatuses.hasAlignedKymo);
            [alignedKymos, alignedKymosStretchFactors, shiftAlignedKymos] = dbmODW.get_aligned_kymos(alignedKymoFileIdxs, alignedKymoFileMoleculeIdxs);
        end

        
        function [moleculeLeftEdgeIdxs, moleculeRightEdgeIdxs] = get_raw_kymo_molecules_edge_idxs(dbmODW, fileIdx, fileMoleculeIdx)
            moleculeStruct = dbmODW.get_molecule_struct(fileIdx, fileMoleculeIdx);
            moleculeLeftEdgeIdxs = moleculeStruct.leftEndArr;
            moleculeRightEdgeIdxs = moleculeStruct.rightEndArr;
        end

        function [] = set_raw_kymo_molecules_edge_idxs(dbmODW, fileIdx, fileMoleculeIdx, moleculeLeftEdgeIdxs, moleculeRightEdgeIdxs)
            moleculeStruct = dbmODW.get_molecule_struct(fileIdx, fileMoleculeIdx);
            moleculeStruct.leftEndArr = moleculeLeftEdgeIdxs;
            moleculeStruct.rightEndArr = moleculeRightEdgeIdxs;
            dbmODW.set_molecule_struct(fileIdx, fileMoleculeIdx, moleculeStruct);
        end
        
        function [moleculesLeftEdgeIdxs, moleculesRightEdgeIdxs] = get_raw_kymos_molecules_edge_idxs(dbmODW, fileIdxs, fileMoleculeIdxs)
            numMolecules = length(fileMoleculeIdxs);
            moleculesLeftEdgeIdxs = cell(numMolecules, 1);
            moleculesRightEdgeIdxs = cell(numMolecules, 1);
            for moleculeNum = 1:numMolecules
                fileIdx = fileIdxs(moleculeNum);
                fileMoleculeIdx = fileMoleculeIdxs(moleculeNum);
                [moleculeLeftEndIdxs, moleculeRightEndIdxs] = dbmODW.get_raw_kymo_molecules_edge_idxs(fileIdx, fileMoleculeIdx);
                moleculesLeftEdgeIdxs{moleculeNum} = moleculeLeftEndIdxs;
                moleculesRightEdgeIdxs{moleculeNum} = moleculeRightEndIdxs;
            end
        end
        
        function [] = set_raw_kymos_molecules_edge_idxs(dbmODW, fileIdxs, fileMoleculeIdxs, moleculesLeftEdgeIdxs, moleculesRightEdgeIdxs)
            numMolecules = length(fileMoleculeIdxs);
            for moleculeNum = 1:numMolecules
                fileIdx = fileIdxs(moleculeNum);
                fileMoleculeIdx = fileMoleculeIdx(moleculeNum);
                moleculeLeftEdgeIdxs = moleculesLeftEdgeIdxs{moleculeNum};
                moleculeRightEdgeIdxs = moleculesRightEdgeIdxs{moleculeNum};
                dbmODW.set_raw_kymo_molecules_edge_idxs(fileIdx, fileMoleculeIdx, moleculeLeftEdgeIdxs, moleculeRightEdgeIdxs);
            end
        end
            
        function [infoScore] = get_info_score(dbmODW, fileIdx, fileMoleculeIdx)
            import OptMap.InfoScore.calc_kymo_info_score;
            import OldDBM.Kymo.Core.find_signal_region_naive;
            infoScore = NaN;
            [moleculeStatus] = dbmODW.get_molecule_statuses(fileIdx, fileMoleculeIdx);
            if (moleculeStatus.hasAlignedKymo && moleculeStatus.passesFilters)
                [alignedKymo, ~, ~] = dbmODW.get_aligned_kymo(fileIdx, fileMoleculeIdx);
                meanAlignedKymoImg = mean(alignedKymo);
                [signalStartIdx, signalEndIdx] = find_signal_region_naive(meanAlignedKymoImg);
                signalRegion = signalStartIdx:signalEndIdx;
                alignedKymoSignalImg = alignedKymo(:, signalRegion);

                % Calculate the information score
                try
                    infoScore = calc_kymo_info_score(alignedKymoSignalImg);
                catch
                end
            end
        end

        function [infoScores] = get_info_scores(dbmODW, fileIdxs, fileMoleculeIdxs)
            numMolecules = length(fileIdxs);
            infoScores = NaN(numMolecules, 1);
            for moleculeNum = 1:numMolecules
                fileIdx = fileIdxs(moleculeNum);
                fileMoleculeIdx = fileMoleculeIdxs(moleculeNum);
                [infoScore] = dbmODW.get_info_score(fileIdx, fileMoleculeIdx);
                infoScores(moleculeNum) = infoScore;
            end
        end
        

        function [moleculeLength] = get_molecule_length(dbmODW, fileIdx, fileMoleculeIdx)
            [moleculeStruct] = dbmODW.get_molecule_struct(fileIdx, fileMoleculeIdx);
            moleculeLength = NaN;
            if isfield(moleculeStruct, 'length') && not(isempty(moleculeStruct.length))
                moleculeLength = moleculeStruct.length;
            end
        end
        
        function [moleculeLengths] = get_molecule_lengths(dbmODW, fileIdxs, fileMoleculeIdxs)
            numMolecules = length(fileIdxs);
            moleculeLengths = NaN(numMolecules, 1);
            for moleculeNum = 1:numMolecules
                fileIdx = fileIdxs(moleculeNum);
                fileMoleculeIdx = fileMoleculeIdxs(moleculeNum);
                [moleculeLength] = dbmODW.get_molecule_length(fileIdx, fileMoleculeIdx);
                moleculeLengths(moleculeNum) = moleculeLength;
            end
        end
        
        
        function [fgKymoTimeAvg] = get_fg_kymo_time_avg(dbmODW, fileIdx, fileMoleculeIdx)
            fgKymoTimeAvg = [];
            [moleculeStruct] = dbmODW.get_molecule_struct(fileIdx, fileMoleculeIdx);
            if isfield(moleculeStruct, 'timeAvKymo')
                fgKymoTimeAvg = moleculeStruct.timeAvKymo;
            end
        end
        
        function [] = set_fg_kymo_time_avg(dbmODW, fileIdx, fileMoleculeIdx, fgKymoTimeAvg)
            moleculeStruct = dbmODW.get_molecule_struct(fileIdx, fileMoleculeIdx);
            moleculeStruct.timeAvKymo = fgKymoTimeAvg;
            dbmODW.set_molecule_struct(fileIdx, fileMoleculeIdx, moleculeStruct); 
        end

        function [fgKymoTimeAvgs] = get_fg_kymo_time_avgs(dbmODW, fileIdxs, fileMoleculeIdxs)
            numMolecules = numel(fileIdxs);
            fgKymoTimeAvgs = cell(numMolecules, 1);
            for moleculeNum=1:numMolecules
                fileIdx = fileIdxs(moleculeNum);
                fileMoleculeIdx = fileMoleculeIdxs(moleculeNum);
                fgKymoTimeAvg = dbmODW.get_fg_kymo_time_avg(fileIdx, fileMoleculeIdx);
                fgKymoTimeAvgs{moleculeNum} = fgKymoTimeAvg;
            end
        end
        
        function [fgKymoTimeAvgs, fgKymoTimeAvgFileIdxs, fgKymoTimeAvgFileMoleculeIdxs] = get_all_existing_fg_kymo_time_avgs(dbmODW)
            [fileIdxs, fileMoleculeIdxs] = dbmODW.get_molecule_idxs();
            [moleculeStatuses] = dbmODW.get_molecule_statuses(fileIdxs, fileMoleculeIdxs);
            fgKymoTimeAvgFileIdxs = fileIdxs(moleculeStatuses.hasFgKymoTimeAvg);
            fgKymoTimeAvgFileMoleculeIdxs = fileMoleculeIdxs(moleculeStatuses.hasFgKymoTimeAvg);
            [fgKymoTimeAvgs] = dbmODW.get_fg_kymo_time_avgs(fgKymoTimeAvgFileIdxs, fgKymoTimeAvgFileMoleculeIdxs);
        end

        function [] = set_fg_kymo_time_avgs(dbmODW, fileIdxs, fileMoleculeIdxs, fgKymoTimeAvgs)
            numFgKymoTimeAvgs = length(fgKymoTimeAvgs);

            for fgKymoTimeAvgNum = 1:numFgKymoTimeAvgs
                fileIdx = fileIdxs(fgKymoTimeAvgNum);
                fileMoleculeIdx = fileMoleculeIdxs(fgKymoTimeAvgNum);
                fgKymoTimeAvg = fgKymoTimeAvgs{fgKymoTimeAvgNum};
                dbmODW.set_fg_kymo_time_avg(fileIdx, fileMoleculeIdx, fgKymoTimeAvg);
            end
        end
    end
    
    methods (Static)
        function [isValid, validationErrMsg] = validate_averaging_window_width(averagingWindowWidth)
            validationErrMsg = '';
            isValid = isscalar(averagingWindowWidth) && ...
                not(...
                    isnan(averagingWindowWidth) || ...
                    (averagingWindowWidth < 1) || ...
                    (rem(averagingWindowWidth, 2) ~= 1));
            if not(isValid)
                validationErrMsg = 'Averaging window width must be a positive odd integer';
            end
        end

        function [averagingWindowWidth] = update_prompt_averaging_window_width(dbmODW)
            % potentially unused?
            
            averagingWindowWidth = prompt_averaging_window_width(dbmODW.get_averaging_window_width());
            dbmODW.set_averaging_window_width(averagingWindowWidth);
        end

        function [averagingWindowWidth] = prompt_averaging_window_width(defaultAveragingWindowWidth)
            import OldDBM.General.DataWrapper;
            
            averagingWindowWidth = defaultAveragingWindowWidth;
            prompt = {'Averaging window width:'};
            defaultVal = {num2str(defaultAveragingWindowWidth)};
            dlg_title = 'Set option';
            num_lines = 1;
            answer = inputdlg(prompt,dlg_title,num_lines,defaultVal);

            if isempty(answer)
                return;
            end

            % Parse the inputs.
            averagingWindowWidth = str2double(answer{1});
            [isValid, validationErrMsg] = DataWrapper.validate_averaging_window_width(averagingWindowWidth);

            if not(isValid)
                disp(validationErrMsg)
                return;
            end
        end
    end
    
    methods (Static, Access = private)

        function [filter_fn] = generate_collective_filter_fn(minMoleculeLengthEstThresh_px, minInfoScoreThresh)
            [mol_length_filter_fn] = generate_mol_length_filter_fn(minMoleculeLengthEstThresh_px);
            [infoscore_filter_fn] = generate_info_score_filter_fn(minInfoScoreThresh);
            function passesFiltersVect = passes_filters(dbmODW, fileIdxs, fileMoleculeIdxs)
                passesFiltersVect = arrayfun(...
                    @(fileIdx, fileMoleculeIdx)...
                        mol_length_filter_fn(dbmODW, fileIdx, fileMoleculeIdx) &...
                        infoscore_filter_fn(dbmODW, fileIdx, fileMoleculeIdx),...
                    fileIdxs, fileMoleculeIdxs ...
                );
            end
            filter_fn = @passes_filters;
            
            function [mol_length_filter_fn] = generate_mol_length_filter_fn(minMoleculeLengthEstThresh_px)
                function passesMoleculeLengthFilter = passes_mol_length_filter(dbmODW, fileIdx, fileMoleculeIdx)
                    moleculeStruct = dbmODW.get_molecule_struct(fileIdx, fileMoleculeIdx);
                    moleculeLengthEstimate_px = moleculeStruct.length;
                    if isempty(moleculeLengthEstimate_px)
                        moleculeLengthEstimate_px = NaN;
                    end
                    passesMoleculeLengthFilter = isnan(moleculeLengthEstimate_px) | (moleculeLengthEstimate_px >= minMoleculeLengthEstThresh_px);

                end
                mol_length_filter_fn = @passes_mol_length_filter;
            end

            function [info_score_filter_fn] = generate_info_score_filter_fn(minInfoScoreThresh)
                function passesInfoScoreFilter = passes_info_score_filter(dbmODW, fileIdx, fileMoleculeIdx)
                    infoScore = dbmODW.get_info_score(fileIdx, fileMoleculeIdx);
                    passesInfoScoreFilter = isnan(infoScore) | (infoScore >= minInfoScoreThresh);
                end
                info_score_filter_fn = @passes_info_score_filter;
            end
        end
    end
    
    methods (Access = private)
    
        function [] = set_averaging_window_width(dbmODW, averagingWindowWidth)
            dbmODW.DBMMainstruct.windowWidth = averagingWindowWidth;
        end
    end
end
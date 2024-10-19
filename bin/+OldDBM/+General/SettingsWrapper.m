classdef SettingsWrapper < handle
    % SETTINGSWRAPPER Functions for maintaining compatibility with
    %  old structure of DBM settings data
    
    properties %(Access = private)
        DBMSettingsstruct
    end
    
    methods
        function [dbmOSW] = SettingsWrapper(settingsStruct)
            dbmOSW.DBMSettingsstruct = settingsStruct;
        end
        
        function [] = update_settings(dbmOSW, dbmOSW2)
            settingsStruct = dbmOSW2.DBMSettingsstruct;
            dbmOSW.DBMSettingsstruct = settingsStruct;
            
            
            fprintf('Loaded the following settings:\n');
            sectionFieldNames = fieldnames(settingsStruct);
            numSections = length(sectionFieldNames);
            for sectionNum = 1:numSections
                sectionFieldName = sectionFieldNames{sectionNum};
                fprintf('[%s]:\n', sectionFieldName);
                disp(settingsStruct.(sectionFieldName));
            end
        end
        
        function [outputDirpath] = get_default_export_dirpath(dbmOSW, dataType)
            switch dataType
                case 'session'
                    outputDirpath = dbmOSW.DBMSettingsstruct.dirs.sessions;
                case 'raw_kymo'
                    outputDirpath = dbmOSW.DBMSettingsstruct.dirs.rawKymos;
                case 'aligned_kymo'
                    outputDirpath = dbmOSW.DBMSettingsstruct.dirs.alignedKymos;
                case 'aligned_kymo_time_avg'
                    outputDirpath = dbmOSW.DBMSettingsstruct.dirs.alignedKymos;
                case 'raw_kymo_center_of_mass'
                    outputDirpath = dbmOSW.DBMSettingsstruct.dirs.fileInfo;
                case 'molecule_analysis'
                     outputDirpath = dbmOSW.DBMSettingsstruct.dirs.analyses;
                case 'raw_kymo_stats'
                     outputDirpath = dbmOSW.DBMSettingsstruct.dirs.stats;
                otherwise
                    import Fancy.AppMgr.AppResourceMgr;
                    appRsrcMgr = AppResourceMgr.get_instance();
                    appDirpath = appRsrcMgr.get_app_dirpath();
                    outputDirpath = appDirpath;
            end
        end
        
        
        function [outputDirpath] = get_default_import_dirpath(dbmOSW, dataType)
            outputDirpath = dbmOSW.get_default_export_dirpath(dataType);
        end
        
        
        function [defaultConsensusDirpath] = get_default_consensus_dirpath(dbmOSW)
            defaultConsensusDirpath = dbmOSW.DBMSettingsstruct.dirs.consensus;
        end
    end

    methods (Static)
        function [defaultSettingsDBM] = generate_default_settings_struct()
            % Default directory path for output files
            %  absolute or relative to directory of DBM_Gui
            defaultSettingsDBM.dirs.outputs = 'OutputFiles';
            
            % Default directory path for specific types of output files
            %  absolute or relative to outputs directory
            defaultSettingsDBM.dirs.sessions = 'Sessions';
            defaultSettingsDBM.dirs.rawKymos = 'RawKymos';
            defaultSettingsDBM.dirs.alignedKymos = 'AlignedKymos';
            defaultSettingsDBM.dirs.analyses = 'AlignedKymos';
            defaultSettingsDBM.dirs.fileInfo = 'FileInfo';
            defaultSettingsDBM.dirs.consensus = 'Consensus';
            defaultSettingsDBM.dirs.stats = 'Stats';
            
            
            
            
            % Whether to prompt the user for
            %  bitmask-related parameters. (1=Yes, 0=No)
            defaultSettingsDBM.bitmasks.promptForBitmaskingParams = 1;
            
            % width of point-spread function (PSF), 
            %  units nm [default: 300 nm]
            defaultSettingsDBM.bitmasks.PSF_width = 300;
            
            % End-cutting parameter. Determines how large "pieces" are
            %  cut out at ends of a time-averaged barcode, units of the
            %  width of the PSF. 
            defaultSettingsDBM.bitmasks.DeltaCut = 3.0;

            % nm per pixel of the CCD camera
            %  [default: 159.2 nm/pixel]
            defaultSettingsDBM.bitmasks.nmPerPixel = 159.2; 
            
            
            
            
            % Whether to prompt the user for a similarity threshold
            %  for barcode clustering every time. (1=Yes, 0=No)
            defaultSettingsDBM.consensus.promptForBarcodeClusterLimit = 1;
            
            % Default similarity threshold for barcode clustering
            defaultSettingsDBM.consensus.barcodeClusterLimit = 0.75;
            
            import OldDBM.General.SettingsWrapper;
            defaultSettingsDBM = SettingsWrapper.process_DBM_settings(defaultSettingsDBM);

        end
        
        function [dbmOSW] = import_dbm_settings_from_ini(dbmSettingsFilepath)
            % IMPORT_DBM_SETTINGS_FROM_INI - reads data from dbm_settings.ini and returns a struct
            %
            % Authors:
            %  Saair Quaderi (refactoring)
            %  Tobias Ambjörnsson
            %
            import OldDBM.General.SettingsWrapper;

            % Read configuration file
            if nargin < 1
                dbmSettingsFilepath = '';
            end
            [settingsStruct, dbmSettingsFilepath] = SettingsWrapper.read_DBM_settings(dbmSettingsFilepath);
            % Remember to change output below if dbm_settings.ini is changed'
            fprintf('Read DBM settings file: ''%s''\n', dbmSettingsFilepath);
        
            dbmOSW = SettingsWrapper(settingsStruct);
        end
        
        function [dbmOSW] = import_dbm_settings_from_session_path(sessionFilepath)
            % IMPORT_DBM_SETTINGS_FROM_SESSION_PATH - reads settings from session path and returns a struct
            %
            % Authors:
            %  Saair Quaderi
            %

            % Read configuration file
            import OldDBM.General.SettingsWrapper;
            
            sessionStruct = load(sessionFilepath, 'DBMSettingsstruct');
            if not(isfield(sessionStruct, 'DBMSettingsstruct'))
                error('Failed to recognize the format of the session file');
            end
            
            settingsStruct = sessionStruct.DBMSettingsstruct;
            dbmOSW = SettingsWrapper(settingsStruct);
        end
        
        
        
        function [defaultSettingsFilepath] = get_default_DBM_ini_filepath()
            defaultSettingsFilename = 'DBM.ini';
            import Fancy.AppMgr.AppResourceMgr;
            defaultSettingsDirpath = AppResourceMgr.get_dirpath('SettingFiles');
            defaultSettingsFilepath = fullfile(defaultSettingsDirpath, defaultSettingsFilename);
        end
        
        function [dbmParamsIniFilepath] = prompt_DBM_ini_filepath(defaultSettingsFilepath)
            import OldDBM.General.SettingsWrapper;
            if nargin < 1
                defaultSettingsFilepath = SettingsWrapper.get_default_DBM_ini_filepath();
            end
            
            [settingsFilename, settingsDirpath] = uigetfile('*.ini', 'Select DBM Settings Ini File', defaultSettingsFilepath);
            if isequal(settingsDirpath, 0)
                dbmParamsIniFilepath = [];
                return;
            end
            dbmParamsIniFilepath = fullfile(settingsDirpath, settingsFilename);
        end
        
        function [fileParamsDBM, dbmParamsIniFilepath] = read_DBM_settings(dbmParamsIniFilepath, appDirpath)
            import OldDBM.General.SettingsWrapper;
            import Fancy.IO.ini2struct;
            
            if (nargin < 1) || isempty(dbmParamsIniFilepath)
                dbmParamsIniFilepath = SettingsWrapper.prompt_DBM_ini_filepath();
            end
            if (nargin < 2) || isempty(appDirpath)
                import Fancy.AppMgr.AppResourceMgr;
                appDirpath = AppResourceMgr.get_app_dirpath();
            end
            fileParamsDBM = SettingsWrapper.process_DBM_settings(ini2struct(dbmParamsIniFilepath), appDirpath);
        end
    end
    methods(Static, Access = private)
        function fileParamsDBM = process_DBM_settings(fileParamsDBM, appDirpath)
            if not(isfield(fileParamsDBM, 'dirs'))
                fileParamsDBM.dirs = [];
            end

            % Default output dir settings

            if not(isfield(fileParamsDBM.dirs, 'outputs'))
                fileParamsDBM.dirs.outputs = 'OutputFiles';
            end
            import Fancy.Utils.resolve_path;
            fileParamsDBM.dirs.outputs = resolve_path(appDirpath, fileParamsDBM.dirs.outputs, true);
            
            import Fancy.IO.mkdirp;
            mkdirp(fileParamsDBM.dirs.outputs)

            outputDirTypes = {'sessions', 'rawKymos' 'alignedKymos' 'analyses', 'fileInfo', 'consensus', 'stats'};
            defaultDefaultOutputDirs = {'Sessions', 'RawKymos', 'AlignedKymos', 'Analyses', 'FileInfo', 'Consensus', 'Stats'};

            for sectionFieldNameIdx=1:length(outputDirTypes)
                outputDirType = outputDirTypes{sectionFieldNameIdx};
                defaultDefaultOutputDir = defaultDefaultOutputDirs{sectionFieldNameIdx};
                if not(isfield(fileParamsDBM.dirs, outputDirType))
                    fileParamsDBM.dirs.(outputDirType) = defaultDefaultOutputDir{sectionFieldNameIdx};
                end
            end


            % Default consensus settings
            if not(isfield(fileParamsDBM, 'consensus'))
                fileParamsDBM.consensus = [];
            end

            if not(isfield(fileParamsDBM.consensus, 'promptForBarcodeClusterLimit'))
                fileParamsDBM.consensus.promptForBarcodeClusterLimit = 0;
            end

            if not(isfield(fileParamsDBM.consensus, 'barcodeClusterLimit'))
                fileParamsDBM.consensus.barcodeClusterLimit = 0.75;
            end

            import Fancy.Utils.resolve_path;
            sectionFieldNames = fieldnames(fileParamsDBM);
            for sectionFieldNameIdx = 1:length(sectionFieldNames) 
                sectionFieldName = sectionFieldNames{sectionFieldNameIdx};
                if (strcmp(sectionFieldName, 'dirs') == 1)
                    % Resolve full paths of directories
                    settingNames = fieldnames(fileParamsDBM.(sectionFieldName));
                    for settingNameIdx=1:numel(settingNames)
                        settingName = settingNames{settingNameIdx};
                        if (any(ismember(outputDirTypes, settingName)))
                            resolvedPath = resolve_path(fileParamsDBM.dirs.outputs, fileParamsDBM.(sectionFieldName).(settingName), true);
                            fileParamsDBM.(sectionFieldName).(settingName) = resolvedPath;
                        end
                    end
                end
            end

            % Create any output directories that may need to be created
            settingNames = fieldnames(fileParamsDBM.dirs);
            for settingNameIdx = 1:length(settingNames)
                settingName = settingNames{settingNameIdx};
                mkdirp(fileParamsDBM.dirs.(settingName))
            end
        end
        
    end
    
end
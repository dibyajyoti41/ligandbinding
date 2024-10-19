function [] = export_dbm_session_struct_mat(dbmODW, dbmOSW, defaultOutputDirpath)
    timestamp = datestr(clock(), 'yyyy-mm-dd_HH_MM_SS');
    defaultOutputFilename = sprintf('session_%s.mat', timestamp);
    defaultOutputFilepath = fullfile(defaultOutputDirpath, defaultOutputFilename);
    [outputMatFilename, outputDirpath] = uiputfile('*.mat', 'Save Session', defaultOutputFilepath);
    if isequal(outputDirpath, 0)
        return;
    end
    outputFilepath = fullfile(outputDirpath, outputMatFilename);
    DBMMainstruct = dbmODW.DBMMainstruct; %#ok<NASGU>
    DBMSettingsstruct = dbmOSW.DBMSettingsstruct; %#ok<NASGU>
    save(outputFilepath, 'DBMMainstruct', 'DBMSettingsstruct');
end
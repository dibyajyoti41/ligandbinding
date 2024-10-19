function [consensusSettings, dbmSettingsStruct] = get_default_consensus_settings(dbmParamsIniFilepath)
    if nargin < 1
        dbmParamsIniFilepath = [];
    end
    
    import OldDBM.General.SettingsWrapper;
    [dbmSettingsStruct] = SettingsWrapper.read_DBM_settings(dbmParamsIniFilepath);
    consensusSettings = dbmSettingsStruct.consensus;
end
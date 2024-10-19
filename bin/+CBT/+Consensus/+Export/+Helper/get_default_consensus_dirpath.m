function defaultConsensusDirpath = get_default_consensus_dirpath(dbmOSW)
    if nargin < 1
        import OldDBM.General.SettingsWrapper;
        defaultSettingsFilepath = SettingsWrapper.get_default_DBM_ini_filepath();
        
        import OldDBM.General.SettingsWrapper;
        dbmOSW = SettingsWrapper.import_dbm_settings_from_ini(defaultSettingsFilepath);
    end
    defaultConsensusDirpath = dbmOSW.get_default_consensus_dirpath();
end
function [aborted, sessionFilepath] = prompt_dbm_session_filepath()
    import Fancy.AppMgr.AppResourceMgr;
    appRsrcMgr = AppResourceMgr.get_instance();
    appDirpath = appRsrcMgr.get_app_dirpath();
    defaultSessionDirpath = appDirpath;
    promptTitle = 'Select a DBM session file for kymo import';

    [sessionFilename, dirpath] = uigetfile({'*.mat;'}, promptTitle, defaultSessionDirpath);
    aborted = isequal(dirpath, 0);
    if aborted
        sessionFilepath = '';
    else
        sessionFilepath = fullfile(dirpath, sessionFilename);
    end
end
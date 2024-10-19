function [sessionFilepath] = try_prompt_single_session_filepath(defaultSessionDirpath)
    if nargin < 1
        import Fancy.AppMgr.AppResourceMgr;
        appRsrcMgr = AppResourceMgr.get_instance();
        appDirpath = appRsrcMgr.get_app_dirpath();
        defaultSessionDirpath = appDirpath;
    end

    [sessionFilename, sessionDirpath] = uigetfile({'*.mat;'}, 'Select session file to import', defaultSessionDirpath);

    if isequal(sessionDirpath, 0)
        sessionFilepath = '';
        return;
    end

    sessionFilepath = fullfile(sessionDirpath, sessionFilename);
end
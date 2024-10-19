function [dbmODW, dbmOSW] = try_loading_from_session_file(sessionFilepath)
    if nargin < 1
        import OldDBM.General.Import.try_prompt_single_session_filepath;
        sessionFilepath = try_prompt_single_session_filepath();

        if isempty(sessionFilepath)
            error('No session file was selected');
        end
    end
    
    hFigsPrev = findall(0, 'Type', 'figure');

    sessionStruct = load(sessionFilepath, 'DBMMainstruct', 'DBMSettingsstruct');
    if (not(isfield(sessionStruct,'DBMMainstruct')) || ...
            not(isfield(sessionStruct, 'DBMSettingsstruct')))
        error('Failed to recognize the format of the session file');
    end

    mainStruct = sessionStruct.DBMMainstruct;
    if isfield(mainStruct, 'fig') && isgraphics(mainStruct.fig, 'figure')
        close(setdiff(mainStruct.fig, hFigsPrev));
    end
    import OldDBM.General.DataWrapper;
    dbmODW = DataWrapper(mainStruct, sessionFilepath);
    

    settingsStruct = sessionStruct.DBMSettingsstruct;
    import OldDBM.General.SettingsWrapper;
    dbmOSW = SettingsWrapper(settingsStruct);
end
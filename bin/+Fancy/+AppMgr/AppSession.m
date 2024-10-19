classdef (Sealed) AppSession < handle
    % APPSESSION - Application session class helps manage
    %   functionality related to an application session
    %
    % Note: An instance of APPSESSION should be retrieved
    %   using get_instance since it will be defined as a persistent
    %   variable using the singleton design pattern.
    %
    % http://mathworks.com/help/matlab/ref/persistent.html
    % http://mathworks.com/help/matlab/matlab_oop/controlling-the-number-of-instances.html
    %
    % Authors:
    %   Saair Quaderi
    
    properties (Constant)
        Version = [0, 0, 1];
    end
        
    properties (Constant, Hidden)
        AppSessionsDataPoolID = sprintf('%sAppSessions', Fancy.AppMgr.AppDataPoolViewScreen.HideDataPoolPrefix);
    end
    properties(SetAccess = private)
        MatlabUser = '';
        EnvUser = '';
        EnvHost = '';
        
        CreationTime = [];
    end
    
    methods (Access = private)
        function [appSess] = AppSession()
            % APPSESSION - Constructor for the APPSESSION object
            %
            % Outputs:
            %   appSess
            %     the application session object
            %
            % Authors:
            %   Saair Quaderi
            
            appSess.CreationTime = clock();
            appSess.update_user_env_info();
        end
    end
    methods
        function [] = update_user_env_info(appSess)
            % UPDATE_USER_ENV_INFO - updates the matlab user, environment
            %   user, and environment host information of the appSess
            %   object
            %
            % Inputs:
            %   appSess
            %     the application session object
            %
            % Authors:
            %   Saair Quaderi
            
            matlabUser = '';
            try
                matlabUser = feval(@(s) s.user, license('inuse', 'MATLAB'));
            catch
            end
            envHostName = '';
            envUser= '';
            ret = 1;
            try
                [ret, envHostName] = system('hostname');   
            catch
            end
            try
                if isunix() 
                    envUser = getenv('USER'); 
                    if ret ~= 0
                        envHostName = getenv('HOSTNAME'); 
                    end
                else 
                    envUser = getenv('username');
                    if ret ~= 0
                        envHostName = getenv('COMPUTERNAME');
                    end
                end
            catch
            end
            matlabUser = strtrim(matlabUser);
            envUser = strtrim(envUser);
            envHostName = strtrim(envHostName);
            
            
            appSess.MatlabUser = matlabUser;
            appSess.EnvUser = envUser;
            appSess.EnvHost = envHostName;
        end
        
        function [sessSummaryPrefixStr, creationTimestamp] = get_session_summary_strs(appSess)
            % GET_SESSION_SUMMARY_STRS - creates the session summary
            %  strings which contains information about the environment,
            %  user, and the time
            %
            % Inputs:
            %   appSess
            %     the application session object
            %
            % Outputs:
            %   sessSummaryPrefixStr
            %     the string concatenating the environment host, environment
            %     user, and matlab user information with underscores
            %   creationTimestamp
            %     the timestamp for the Creation Time of the application
            %     session
            %
            % Authors:
            %   Saair Quaderi
            
            envHost = appSess.EnvHost;
            envUser = appSess.EnvUser;
            matlabUser = appSess.MatlabUser;
            sessSummaryPrefixStr = sprintf('%s_%s_%s', envHost, envUser, matlabUser);
            
            creationTimestamp = '';
            if not(isempty(appSess.CreationTime))
                creationTimestamp = datestr(appSess.CreationTime, 'yyyy-mm-dd HH:MM:SS');
            end
        end
    end
    methods (Static)
        function appSess = get_instance(varargin)
            % GET_INSTANCE - returns an instance of AppSession,
            %   generating it if a valid one doesn't already exist and
            %   returning the valid existing one otherwise
            %
            % Inputs:
            %  varargin (any number of arguments)
            %   If one of arguments is 'verbose' (case insensitive), the
            %    function prints vebose console text about what occurs
            %    (e.g. specifying whether the app session was
            %      retrieved or newly generated)
            %
            % Outputs:
            %  appSess
            %    the instance of AppSession
            %
            % Authors:
            %   Saair Quaderi
            
            strArgMask = cellfun(@(arg) (ischar(arg) && isrow(arg)), varargin);
            strArgs = cellfun(@lower, varargin(strArgMask), 'UniformOutput', false);
            verbose = not(isempty(intersect(strArgs, {'verbose'})));
            
            persistent persistentLocalAppSess
            import Fancy.AppMgr.AppSession;
            if isempty(persistentLocalAppSess) || not(isvalid(persistentLocalAppSess))
                
                appSess = AppSession();
                import Fancy.AppMgr.AppSession;
                if not(AppDataPoolMgr.has_data_pool(AppSession.AppSessionsDataPoolID))
                    AppDataPoolMgr.create_new_data_pool(AppSession.AppSessionsDataPoolID);
                end

                [sessSummaryPrefixStr, timestamp] = appSess.get_session_summary_strs();
                sessDataID = sprintf('%s (%s)', sessSummaryPrefixStr, timestamp);

                import Fancy.AppMgr.AppDataPoolMgr;
                appDataPoolMgr = AppDataPoolMgr.get_instance();
                [appSessions, appSessionItemIDs] = appDataPoolMgr.get_data_items(AppSession.AppSessionsDataPoolID);

                relatedSessionMask = strncmp(sessSummaryPrefixStr, appSessionItemIDs(:), length(sessSummaryPrefixStr));
                numRelatedSessions = sum(relatedSessionMask);
                if (numRelatedSessions == 0)
                    appDataPoolMgr.update_data_item(AppSession.AppSessionsDataPoolID, sessDataID, appSess);
                    persistentLocalAppSess = appSess;

                    if verbose
                        fprintf('Using Newly Generated App Session: %s\n', sessDataID);
                    end
                    return;
                end
                relatedSessionIDs = appSessionItemIDs(relatedSessionMask);
                relatedSessions = appSessions(relatedSessionMask);
                
                import Fancy.UI.FancyInput.dropdown_dialog;
                [~, idx] = dropdown_dialog('App Session Selection', 'Choose a session to use', relatedSessionIDs);
                appSess = relatedSessions{idx};
                persistentLocalAppSess = appSess;
            else
                appSess = persistentLocalAppSess;
                if verbose
                    [sessSummaryPrefixStr, timestamp] = appSess.get_session_summary_strs();
                    sessDataID = sprintf('%s (%s)', sessSummaryPrefixStr, timestamp);
                    fprintf('Retrieved App Session from Memory: %s\n', sessDataID);
                end
            end
        end
    end
end

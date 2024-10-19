classdef (Sealed) AppResourceMgr < handle
    % APPRESOURCEMGR - Application resource management class helps manage
    %   functionality related to resources available to the app
    %
    % Note: An instance of APPRESOURCEMGR should be retrieved
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
        
        AppResourceMgrsDataPoolID = sprintf('%sAppResourceMgrs', Fancy.AppMgr.AppDataPoolViewScreen.HideDataPoolPrefix);
    end
    properties (Constant, Access = private)
        RequiredMatlabDirname = 'MATLAB'
        RequiredMatlabParentDirname = 'src'
        
        IVY_GROUP_ID = 'org.apache.ivy';
        IVY_ARTIFACT_ID = 'ivy';
        IVY_VERSION = '2.4.0';
        IVY_EXPECTED_CHECKSUM = '8C88B943FCD643D5E592B86179C6FBEB'; % MD5 Checksum (Note that this is not resistant to collision attacks)
        IVY_CLASS = 'Ivy';
    end
    properties
        AppSession
        Dirpaths = containers.Map('KeyType', 'char', 'ValueType', 'char');
    end
    
    methods (Access = private)
        function [appRsrcMgr] = AppResourceMgr(appSess)
            % APPRESOURCEMGR - Constructor for instances of APPRESOURCEMGR
            %
            %  Inputs:
            %    appSess (optional; defaults to instance of AppSession
            %     returned by AppSession.get_instance)
            %
            %  Outputs:
            %    appRsrcMgr
            %      AppResourceMgr object
            % 
            %  Authors:
            %    Saair Quaderi
            
            if nargin < 1
                appSess = [];
            end
            
            if isempty(appSess)
                import Fancy.AppMgr.AppSession;
                appSess = AppSession.get_instance();
            end
            appRsrcMgr.AppSession = appSess;
            
            
            import Fancy.AppMgr.AppResourceMgr;
            
            matlabPackagesDirpath = AppResourceMgr.find_root_package_parent_dirpath();
            [srcDirpath, matlabPackagesDirname] = fileparts(matlabPackagesDirpath);
            [appDirpath, srcDirname] = fileparts(srcDirpath);
            if not(strcmp(matlabPackagesDirname,  AppResourceMgr.RequiredMatlabDirname))
                error('Expected matlab packages directory to be named ''%s''', ...
                    AppResourceMgr.RequiredMatlabDirname);
            end
            if not(strcmp(srcDirname,  AppResourceMgr.RequiredMatlabParentDirname))
                error('Expected matlab packages directory''s parent directory to be named ''%s''', ...
                    AppResourceMgr.RequiredMatlabParentDirname);
            end
            
            % Have to call the non-static private p_* versions because
            % calling the static versions could retrigger AppResourceMgr
            % constructor through get_app_instance which might result in an
            % infinite loop
            appRsrcMgr.p_init_dirpath('App', appDirpath);
            appRsrcMgr.p_init_dirpath('bin', fullfile(appDirpath, 'bin'));
            appRsrcMgr.p_init_dirpath('lib/jar', fullfile(appDirpath, 'lib', 'jar'));
            appRsrcMgr.p_init_dirpath('MATLAB', matlabPackagesDirpath);
            appRsrcMgr.p_init_dirpath('DataCache', fullfile(appDirpath, 'DataCache'));
            appRsrcMgr.p_init_dirpath('SettingFiles', fullfile(appDirpath, 'SettingFiles'));
            
            appRsrcMgr.p_add_all_packages();
            appRsrcMgr.p_download_ivy_jar();
        end
        
        function [] = p_add_all_packages(appRsrcMgr)
            % P_ADD_ALL_PACKAGES - Private function to add the package
            %   directory path and all descendant subdirectory paths to the
            %   Matlab path
            %
            % Inputs:
            %   appRsrcMgr
            %     the Application Resource Manager object which contains
            %     the information about where the packages directory is
            %     located
            %
            % Side effects:
            %    the Matlab path has the package directory path and all
            %    descendant subdirectory paths added to it
            %
            % Authors:
            %   Saair Quaderi
            
            addpath(genpath(appRsrcMgr.p_get_matlab_dirpath()));
        end
        
        function [] = p_init_dirpath(appRsrcMgr, dirpathKey, dirpath)
            % P_INIT_DIRPATH - Private function to  making the directory 
            %   if it does not already exist avd add/set the directory
            %   path specified in the Application Resource Manager with the
            %   specified path key
            %
            % Inputs:
            %   appRsrcMgr
            %     the Application Resource Manager object which contains
            %     the information about where the packages directory is
            %     located
            %   dirpathKey
            %     the key to associate with the directory path
            %   dirpath
            %     the directory path
            %
            % Authors:
            %   Saair Quaderi
            
            import Fancy.IO.mkdirp;
            mkdirp(dirpath);
            if not(isKey(appRsrcMgr.Dirpaths, dirpathKey))
                appRsrcMgr.p_set_dirpath(dirpathKey, dirpath);
            else
                error('Path already exists');
            end
        end
        
        function [] = p_set_dirpath(appRsrcMgr, dirpathKey, dirpath)
            % P_SET_DIRPATH - Private function to set the directory
            %   path specified in the Application Resource Manager with the
            %   specified path key
            %
            % Inputs:
            %   appRsrcMgr
            %     the Application Resource Manager object which contains
            %     the information about where the packages directory is
            %     located
            %   dirpathKey
            %     the key to associate with the directory path
            %   dirpath
            %     the directory path
            %
            % Authors:
            %   Saair Quaderi
            
            appRsrcMgr.Dirpaths(dirpathKey) = dirpath;
        end
        
        function [appDirpath] = p_get_app_dirpath(appRsrcMgr)
            % P_GET_APP_DIRPATH - Private function to get the directory
            %   path for the application
            %
            % Inputs:
            %   appRsrcMgr
            %     the Application Resource Manager object which contains
            %     the information about where the application directory is
            %     located
            %
            % Outputs:
            %   appDirpath
            %      the directory path for the application
            %
            % Authors:
            %   Saair Quaderi
            
            appDirpath = appRsrcMgr.Dirpaths('App');
        end
        
        function [matlabDirpath] = p_get_matlab_dirpath(appRsrcMgr)
            % P_GET_PACKAGES_DIRPATH - Private function to get the
            %   directory path for MATLAB code
            %
            % Inputs:
            %   appRsrcMgr
            %     the Application Resource Manager object which contains
            %     the information about where the packages directory is
            %     located
            %
            % Outputs:
            %   matlabDirpath
            %      the directory path for the packages
            %
            % Authors:
            %   Saair Quaderi
            
            matlabDirpath = appRsrcMgr.Dirpaths('MATLAB');
        end
        
        function [libJarDirpath] = p_get_lib_jar_dirpath(appRsrcMgr)
            % P_GET_LIB_JAR_DIRPATH - Private function to get the
            %   directory path for library jar files
            %
            % Inputs:
            %   appRsrcMgr
            %     the Application Resource Manager object which contains
            %     the information about where the lib/jar/ directory is
            %     located
            %
            % Outputs:
            %   libJarDirpath
            %      the directory path for the library jar files
            %
            % Authors:
            %   Saair Quaderi
            
            libJarDirpath = appRsrcMgr.Dirpaths('lib/jar');
        end
        
        function [dirpath] = p_get_dirpath(appRsrcMgr, dirpathKey)
            % P_GET_DIRPATH - Private function to get the
            %   path associated with the key specified
            %
            % Inputs:
            %   appRsrcMgr
            %     the Application Resource Manager object which contains
            %     the information about where the packages directory is
            %     located
            %   dirpathKey
            %     the key for the directory path of interest
            %
            % Outputs:
            %   dirpath
            %      the directory path for the key of interest
            %
            % Authors:
            %   Saair Quaderi
            
            appDirpath = appRsrcMgr.p_get_app_dirpath();
            if appRsrcMgr.Dirpaths.isKey(dirpathKey)
                dirpath = appRsrcMgr.Dirpaths(dirpathKey);
            else
                defaultDirpath = fullfile(appDirpath, dirpathKey);
                if not(exist(defaultDirpath, 'dir'))
                    defaultDirpath = appDirpath;
                end
                dirpath = uigetdir(defaultDirpath, sprintf('Select ''%s'' directory', dirpathKey));
                if not(isequal(dirpath, 0))    
                    appRsrcMgr.Dirpaths(dirpathKey) = dirpath;
                end
            end
        end
        
        
        function [] = p_add_lib_jar_paths(appRsrcMgr)
            % P_ADD_LIB_JAR_PATHS - Private function to add the
            %   jar files directly in the app's lib/jar/ directory
            %   to the 
            %
            % Inputs:
            %   appRsrcMgr
            %     the Application Resource Manager object which contains
            %     the information about where the app directory is
            %     located
            %
            % Side Effects:
            %   dynamic javaclasspath is updated
            %
            % Authors:
            %   Saair Quaderi
            
            libJarDirpath = appRsrcMgr.p_get_lib_jar_dirpath();
            dirStruct = dir(libJarDirpath);
            dirFilenames = feval(@(x) x(:), {dirStruct(~[dirStruct.isdir]).name});
            [~, ~, dirFileExts] = cellfun(@fileparts, dirFilenames, 'UniformOutput', false);
            dirJarFilepaths = fullfile(libJarDirpath, dirFilenames(cellfun(@(ext) strcmpi(ext, '.jar'), dirFileExts)));
            currDynamicJavaPaths = javaclasspath('-dynamic');
            newDynamicJarFilepaths = setdiff(dirJarFilepaths, currDynamicJavaPaths);
            if not(isempty(newDynamicJarFilepaths))
                javaaddpath(newDynamicJarFilepaths);
            end
        end
        
        function [] = p_download_ivy_jar(appRsrcMgr, verifyIvyCheckSumTF)
            % P_DOWNLOAD_IVY_JAR - Private function to download the
            %   jar file for Apache Ivy (which can then be used to
            %   facilitate downloads of other jar dependencies)
            %
            % Inputs:
            %   appRsrcMgr
            %     the Application Resource Manager object which contains
            %     the information about where the jar directory is
            %     located
            %   verifyIvyCheckSumTF
            %     true if the download's MD5 checksum should be checked
            %      once the file is downloaded
            %
            % Side Effects:
            %   Downloads a jar file into the jar directory,
            %    assuming the download goes smoothly 
            %
            % Authors:
            %   Saair Quaderi
            if nargin < 2
                verifyIvyCheckSumTF = true;
            end
            
            libJarDirpath = appRsrcMgr.p_get_lib_jar_dirpath();

            ivyJarFilename = sprintf('%s-%s.jar', AppResourceMgr.IVY_ARTIFACT_ID, AppResourceMgr.IVY_VERSION);
            ivyJarFilepath = fullfile(libJarDirpath, ivyJarFilename);
            import Fancy.AppMgr.AppResourceMgr;
            if not(exist(ivyJarFilepath, 'file'))
                mavenJarBaseURLDirpath = 'http://central.maven.org/maven2/';
                ivyURLDirpath = sprintf('%s%s/%s/%s/', mavenJarBaseURLDirpath, strrep(AppResourceMgr.IVY_GROUP_ID, '.', '/'), AppResourceMgr.IVY_ARTIFACT_ID, AppResourceMgr.IVY_VERSION);
                ivyURL = [ivyURLDirpath, ivyJarFilename];

                try
                    fprintf('Downloading Ivy (''%s'')\n  to ''%s''\n  from ''%s''...\n', ivyJarFilename, libJarDirpath, ivyURLDirpath);
                    websave(ivyJarFilepath, ivyURL);
                catch
                end
            else
                return;
                % % Commenting since this is rather slow and probably
                % unnecessa
                % if verifyIvyCheckSumTF
                %     checksum = Simulink.getFileChecksum(ivyJarFilepath);
                %     if strcmp(checksum, AppResourceMgr.IVY_EXPECTED_CHECKSUM)
                %         fprintf('Ivy jar checksum validated\n');
                %         return;
                %     end
                % end
            end
            
            if exist(ivyJarFilepath, 'file')
                if verifyIvyCheckSumTF
                    checksum = Simulink.getFileChecksum(ivyJarFilepath);
                    if strcmp(checksum, AppResourceMgr.IVY_EXPECTED_CHECKSUM)
                        fprintf('Ivy download attempt completed successfully\n');
                    else
                        fprintf('Ivy jar download had a bad checksum, deleting file\n');
                        delete(ivyJarFilepath);
                    end
                else
                    fprintf('Download attempt completed\n');
                end
            else
                fprintf('Download attempt failed\n');
            end
        end
        
        function [] = p_download_java_dependency(appRsrcMgr, groupId, artifactId, version)
            % P_DOWNLOAD_JAVA_DEPENDENCY - Private function to download the
            %   jar files associated with the specified item
            %   to the app's lib/jar/ directory using Ivy
            %
            % Inputs:
            %   appRsrcMgr
            %     the Application Resource Manager object which contains
            %     the information about where the app directory is
            %     located
            %   groupId
            %     Maven groupId that specifies a project package 
            %   artifactId
            %     Maven artifactId that specifies the name of the jar without
            %      version 
            %   version
            %     Maven version for the version of the jar file
            %
            % Side Effects:
            %   Downloads jar files into the jar directory,
            %    assuming the downloads go smoothly 
            %
            % Authors:
            %   Saair Quaderi
            
            import Fancy.AppMgr.AppResourceMgr;
            ivyFullClassPath = sprintf('%s.%s', AppResourceMgr.IVY_GROUP_ID, AppResourceMgr.IVY_CLASS);
            libJarDirpath = appRsrcMgr.p_get_lib_jar_dirpath();
            ivyJarFilename = sprintf('%s-%s.jar', AppResourceMgr.IVY_ARTIFACT_ID, AppResourceMgr.IVY_VERSION);
            
            appRsrcMgr.add_lib_jar_paths();
            if not(exist(ivyFullClassPath, 'class'))
                appRsrcMgr.download_ivy_jar(false);
                appRsrcMgr.add_lib_jar_paths();
            end
            if not(exist(ivyFullClassPath, 'class'))
                error('''%s'' is missing and is necessary to download dependencies', ivyFullClassPath);
            else
                system(retrieve_ivy_system_cmd(libJarDirpath, ivyJarFilename, groupId, artifactId, version));
                appRsrcMgr.add_lib_jar_paths();
            end
            function sysCmdTxt = retrieve_ivy_system_cmd(libJarDirpath, ivyJarFilename, groupId, artifactId, version)
                ivyJarFilepath = fullfile(libJarDirpath, ivyJarFilename);
                dependencyStr = sprintf('%s %s %s', groupId, artifactId, version);
                retrieveStr = sprintf('"%s[artifact]-[revision](-[classifier]).[ext]"', [libJarDirpath, filesep]);
                sysCmdTxt = sprintf('java -jar %s -dependency %s -retrieve %s', ivyJarFilepath, dependencyStr, retrieveStr);
            end
        end
    end
    methods (Static)
        function [] = init_dirpath(dirpathKey, dirpath)
            % INIT_DIRPATH - Public static function to  making the
            %   directory if it does not already exist avd add/set the
            %   directory path specified in the instance of AppResourceMgr
            %   returned by AppResourceMgr.get_instance
            %
            % Inputs:
            %   dirpathKey
            %     the key to associate with the directory path
            %   dirpath
            %     the directory path
            %
            % Authors:
            %   Saair Quaderi
            
            import Fancy.AppMgr.AppResourceMgr;
            appRsrcMgr = AppResourceMgr.get_instance();
            appRsrcMgr.p_init_dirpath(dirpathKey, dirpath);
        end
        
        function [] = set_dirpath(dirpathKey, dirpath)
            % SET_DIRPATH - Public static function to set the directory
            %   path specified in the instance of AppResourceMgr returned
            %   by AppResourceMgr.get_instance
            %
            % Inputs:
            %   dirpathKey
            %     the key to associate with the directory path
            %   dirpath
            %     the directory path
            %
            % Authors:
            %   Saair Quaderi
            
            import Fancy.AppMgr.AppResourceMgr;
            appRsrcMgr = AppResourceMgr.get_instance();
            appRsrcMgr.p_set_dirpath(dirpathKey, dirpath);
        end
        
        function [dirpath] = get_dirpath(dirpathKey)
            % GET_DIRPATH - Public static function to get the
            %   path associated with the key specified from the instance of
            %   AppResourceMgr returned by AppResourceMgr.get_instance
            %
            % Inputs:
            %   dirpathKey
            %     the key for the directory path of interest
            %
            % Outputs:
            %   dirpath
            %      the directory path for the key of interest
            %
            % Authors:
            %   Saair Quaderi
            import Fancy.AppMgr.AppResourceMgr;
            appRsrcMgr = AppResourceMgr.get_instance();
            dirpath = appRsrcMgr.p_get_dirpath(dirpathKey);
        end

        function [matlabDirpath] = get_matlab_dirpath()
            % GET_MATLAB_DIRPATH - Public static function to get the
            %   directory path for MATLAB packages from the instance of
            %   AppResourceMgr returned by AppResourceMgr.get_instance
            %
            % Outputs:
            %   matlabDirpath
            %      the directory path for the MATLAB code
            %
            % Authors:
            %   Saair Quaderi
            
            import Fancy.AppMgr.AppResourceMgr;
            appRsrcMgr = AppResourceMgr.get_instance();
            matlabDirpath = appRsrcMgr.p_get_matlab_dirpath();
        end
        
        function [appDirpath] = get_app_dirpath()
            % GET_APP_DIRPATH - Public static function to get the
            %   directory path for the application from the instance of
            %   AppResourceMgr returned by AppResourceMgr.get_instance
            %
            % Outputs:
            %   appDirpath
            %      the directory path for the application
            %
            % Authors:
            %   Saair Quaderi
            
            import Fancy.AppMgr.AppResourceMgr;
            appRsrcMgr = AppResourceMgr.get_instance();
            appDirpath = appRsrcMgr.p_get_app_dirpath();
        end
        
        function [] = add_all_packages()
            % ADD_ALL_PACKAGES - Public static function to add the package
            %   directory path and all descendant subdirectory paths to the
            %   Matlab path from the instance of AppResourceMgr returned by
            %   AppResourceMgr.get_instance
            %
            % Side effects:
            %    the Matlab path has the package directory path and all
            %    descendant subdirectory paths added to it
            %
            % Authors:
            %   Saair Quaderi
            
            import Fancy.AppMgr.AppResourceMgr;
            appRsrcMgr = AppResourceMgr.get_instance();
            appRsrcMgr.p_add_all_packages();
        end
        
        function [] = add_lib_jar_paths()
            % ADD_LIB_JAR_PATHS - Private function to add the
            %   jar files directly in the app's lib/jar/ directory
            %   to the 
            %
            % Inputs:
            %   appRsrcMgr
            %     the Application Resource Manager object which contains
            %     the information about where the app directory is
            %     located
            %
            % Side Effects:
            %   javaclasspath 
            %
            % Authors:
            %   Saair Quaderi
            
            import Fancy.AppMgr.AppResourceMgr;
            appRsrcMgr = AppResourceMgr.get_instance();
            appRsrcMgr.p_add_lib_jar_paths();
        end
        
        function appRsrcMgr = get_instance(varargin)
            % GET_INSTANCE - returns an instance of AppResourceMgr,
            %   generating it if a valid one doesn't already exist and
            %   returning the valid existing one otherwise
            %
            % Inputs:
            %  varargin (any number of arguments)
            %   If one of arguments is 'verbose' (case insensitive), the
            %    function prints vebose console text about what occurs
            %    (e.g. specifying whether the app resource manager was
            %      retrieved or newly generated)
            %
            % Outputs:
            %  appRsrcMgr
            %    the instance of AppResourceMgr
            %
            % Authors:
            %   Saair Quaderi
            
            import Fancy.AppMgr.AppSession;
            appSess = AppSession.get_instance();
            
            strArgMask = cellfun(@(arg) (ischar(arg) && isrow(arg)), varargin);
            strArgs = cellfun(@lower, varargin(strArgMask), 'UniformOutput', false);
            verbose = not(isempty(intersect(strArgs, {'verbose'})));
            persistent persistentLocalAppRsrcMgr
            import Fancy.AppMgr.AppResourceMgr;
            if isempty(persistentLocalAppRsrcMgr) || not(isvalid(persistentLocalAppRsrcMgr))
                if verbose
                    fprintf('Generating New App Resource Manager...\n');
                end
                persistentLocalAppRsrcMgr = AppResourceMgr(appSess);
            else
                if verbose
                    fprintf('Retrieving App Resource Manager from Memory...\n');
                end
            end
            appRsrcMgr = persistentLocalAppRsrcMgr;
        end
        
        function [] = download_ivy_jar(verifyIvyCheckSumTF)
            % DOWNLOAD_IVY_JAR - Public static function to download the
            %   jar file for Apache Ivy (which can then be used to
            %   facilitate downloads of other jar dependencies)
            %
            % Inputs:
            %   verifyIvyCheckSumTF (optional; defaults to true)
            %     true if the download's MD5 checksum should be checked
            %      once the file is downloaded
            %
            % Side Effects:
            %   Downloads a jar file into the jar directory,
            %    assuming the download goes smoothly 
            %
            % Authors:
            %   Saair Quaderi
            if nargin < 1
                verifyIvyCheckSumTF = true;
            end
            
            import Fancy.AppMgr.AppResourceMgr;
            appRsrcMgr = AppResourceMgr.get_instance();
            appRsrcMgr.p_download_ivy_jar(verifyIvyCheckSumTF);
        end

        function [] = download_java_dependency(groupId, artifactId, version)
            % DOWNLOAD_JAVA_DEPENDENCY - Public static function to download the
            %   jar files associated with the specified item
            %   to the app's lib/jar/ directory using Ivy
            %
            % Inputs:
            %   groupId
            %     Maven groupId that specifies a project package 
            %   artifactId
            %     Maven artifactId that specifies the name of the jar without
            %      version 
            %   version
            %     Maven version for the version of the jar file
            %
            % Side Effects:
            %   Downloads jar files into the jar directory,
            %    assuming the downloads go smoothly 
            %
            % Authors:
            %   Saair Quaderi

            import Fancy.AppMgr.AppResourceMgr;
            appRsrcMgr = AppResourceMgr.get_instance();
            appRsrcMgr.p_download_java_dependency(groupId, artifactId, version);
        end
    end
    
    methods (Static, Access = private)
        function rootPkgParentDirpath = find_root_package_parent_dirpath()
            % FIND_ROOT_PACKAGE_PARENT_DIRPATH - Private static function
            %  that returns the path of the parent directory of the
            %  root package directory containing this file
            %
            % Outputs:
            %   rootPkgParentDirpath
            %     the path of the parent directory of the root package
            %     directory containing this file
            %
            % Authors:
            %   Saair Quaderi
            
            currFilepath = mfilename('fullpath');
            [currDirpath, ~, ~] = fileparts(currFilepath);
            splitPathStrs = strsplit(fileparts(currDirpath), filesep);
            packageDirStrIdx = length(splitPathStrs);
            while ((packageDirStrIdx > 1) && ...
                    not(isempty(splitPathStrs{packageDirStrIdx})) && ...
                    (strcmp(splitPathStrs{packageDirStrIdx}(1), '+')))
                packageDirStrIdx = packageDirStrIdx - 1;
            end
            rootPkgParentDirpath = fullfile(splitPathStrs{1:packageDirStrIdx});
            if (isunix() && not(isempty(rootPkgParentDirpath)))
                rootPkgParentDirpath = [filesep(), rootPkgParentDirpath];
            end
        end
    end
end

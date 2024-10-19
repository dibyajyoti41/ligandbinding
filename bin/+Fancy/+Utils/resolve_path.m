function resolvedPath = resolve_path(workingDir, specifiedPath, isDirPath)
    if not(isempty(workingDir))
        validateattributes(workingDir, {'char'}, {'vector'}, 1);
    end

    pathLen = length(specifiedPath);
    if pathLen == 0
        resolvedPath = workingDir;
        return
    else
        validateattributes(specifiedPath, {'char'}, {'vector'}, 2);
    end

    if nargin < 3
        isDirPath = false;
    else
        validateattributes(isDirPath, {'logical'}, {'scalar'}, 3);
    end

    resolvedPath = specifiedPath;
    pcSlash = '\';
    unixSlash = '/';
    if ispc
        if (length(strfind(resolvedPath(1:min(pathLen, 3)),':\')) == 1)
            return; %absolute path provided
        end
        resolvedPath = strrep(resolvedPath, unixSlash, filesep); %slashes should be windows-style
    else
        if strcmp(resolvedPath(1), unixSlash)
            return; %absolute path provided
        end
        resolvedPath = strrep(resolvedPath, pcSlash, filesep); %slashes should be unix-style
    end
    resolvedPath = fullfile(workingDir, resolvedPath);
    if (isDirPath && not(isempty(resolvedPath)) && not(resolvedPath(end) == filesep))
        resolvedPath = [resolvedPath, filesep];
    end
end
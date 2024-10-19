function [] = mkdirp(dirpath, msg)
    % MKDIRP - creates directories along a directory path if they
    %     do not exist (and potentially displays a message if it does so)  
    % Inputs:
    %   dirpath
    %     the path to the directory to be created
    %   msg (optional)
    %     message to display when creating a directory
    %
    % Side-effects:
    %   Creates certain directories if they do not exist
    %
    % Authors:
    %   Saair Quaderi
    
    if (exist(dirpath, 'dir') == 0)
        if (nargin < 2)
            msg = ['Creating directory path: ''', dirpath, ''''];
        end
        if not(isempty(msg))
            disp(msg);
        end
        mkdir(dirpath);
    end
end
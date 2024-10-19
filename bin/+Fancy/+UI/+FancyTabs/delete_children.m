function [] = delete_children(hParent)
    % DELETE_CHILDREN - deletes the children of provided parentz
    %
    % Inputs:
    %  hParent:
    %   the handle of the parent whose children should be deleted
    %
    % Side-effects:
    %   deletes the children of provided parents
    %
    % Authors:
    %   Saair Quaderi
    
    hChildren = get(hParent, 'Children');
    numChildren = numel(hChildren);
    for childNum=numChildren:-1:1
        delete(hChildren(childNum));
    end
end
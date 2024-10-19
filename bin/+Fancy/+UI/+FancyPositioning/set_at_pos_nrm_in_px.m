function [hGO, posPx] = set_at_pos_nrm_in_px(hGO, relPosNrm, fnPosPxPostprocess)
    % SET_AT_POS_NRM_IN_PX - sets a graphics object's position in
    %   it's parent in normalized units (converted to pixels) and
    %   then potentially adjusted in pixels
    %
    % Inputs:
    %  hGO
    %    the handle of the graphics object
    %  relPosNrm (optional & skippable)
    %    matlab position in normalized units (pre-px adjustment) where
    %      the object should be relative to its parent
    %      specified in the usual [left bottom width height]
    %      defaulted to [0 0 1 1]
    %  fnPosPxPostprocess (optional)
    %    function that takes the would-be pixelposition of the
    %    graphics object positioned in normalized units if position
    %    was set with relPosNrm and returns the pixel-wise adjusted
    %    position for the graphics object also in pixel units (note
    %    that the height and width will be given a value of 1 if
    %    they are less than 1)
    %
    % Outputs:
    %  posPx
    %    the position in pixels for the graphics object
    %
    %  Side-effects
    %    sets the graphics object at the pixelposition indicated by
    %    posPx. Note that the units (normalized/pixels/etc.) are
    %    converted back such that it remains the same for the
    %    graphics object after this function as it did before it
    %
    % Authors:
    %   Saair Quaderi
    
    import Fancy.UI.FancyPositioning.set_position;
    defaultRelPosNrm = [0, 0, 1, 1];
    if nargin < 3
        fnPosPxPostprocess = @(x) x;
    end
    if (nargin == 2) && isa(relPosNrm, 'function_handle') % 2nd param was skipped
        fnPosPxPostprocess = relPosNrm;
        relPosNrm = defaultRelPosNrm;
    elseif (nargin < 2) || isempty(relPosNrm)
        relPosNrm = defaultRelPosNrm;
    end

    parentPosPx = getpixelposition(get(hGO, 'parent'));
    posPx = (relPosNrm .* [parentPosPx(3:4), parentPosPx(3:4)]);
    posPx = fnPosPxPostprocess(posPx);
    posPx(3:4) = max(1, posPx(3:4)); %minimal height/width of 1px

    set_position(hGO, posPx);
end
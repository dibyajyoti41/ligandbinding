function [] = set_position(hGO, position, tmpUnits, outerposition)
    % SET_POSITION - sets a graphics object to a position (with a
    %    given set of units or defaulting to 'pixels') without
    %    making any lasting change to the object's units
    %
    % Inputs:
    %  hGO
    %    a graphics object handle
    %
    %  position
    %    the position to set the object at specified in the format
    %     [left, bottom, width, height] (the matlab standard)
    %
    %  tmpUnits (optional)
    %    the units with which to set the position
    %     defaults to 'pixels'
    %
    %  outerposition (optional)
    %    if true, sets outerposition instead of position
    %     defaults to false
    %
    % Side-effects:
    %   Position change for graphics object
    %
    % Authors:
    %   Saair Quaderi

    if not(isequal(isgraphics(hGO), true))
        error('Input must be a scalar graphic object handle');
    end

    if (nargin < 3)
        tmpUnits = 'pixels';
        setpixelposition(hGO, position);
    end

    if (nargin < 4)
        outerposition = false;
    end
    if outerposition
        propertyName = 'OuterPosition';
    else
        propertyName = 'Position';
    end

    oldUnits = get(hGO, 'Units');
    if strcmpi(oldUnits, tmpUnits)
        set(hGO, propertyName, position);
    else
        set(hGO, 'Units', tmpUnits, propertyName, position);
        set(hGO, 'Units', oldUnits);
    end
end

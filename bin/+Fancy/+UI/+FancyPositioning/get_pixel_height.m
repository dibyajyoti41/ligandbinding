function heightPx = get_pixel_height(goHandle)
    % GET_PIXEL_HEIGHT - returns the height of the graphics object in
    %   pixels
    %
    % Inputs:
    %   goHandle
    %     the handle of the graphics object
    %
    % Outputs:
    %   heightPx
    %     the height of the graphics object in pixels
    %
    % Authors:
    %   Saair Quaderi
    posPx = getpixelposition(goHandle);
    heightPx = posPx(4);
end

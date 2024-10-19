function [blurredImg] = apply_gaussian_blur(img, squareSmoothingWindowLen_pixels, blurSigmaWidth_pixels)
    % APPLY_GAUSSIAN_BLUR - smooths imgArr using the windowsize 
    %	smoothWindow (see the WPAlign paper for details)
    %
    % Inputs: 
    %	img
    %     the input image to be blurred
    %	squareSmoothingWindowLen_pixels
    %     the length of the sides of the square smoothing window
    %     (in pixels)
    %   blurSigmaWidth_pixels
    %     the length of a single standard deviation of the normal
    %     distribution with which the image will be blurred
    %    (in pixels)
    %
    % Outputs:
    %	blurredImg
    %     the image after gaussian blurring
    % 
    % Authors:
    %   Charleston Noble
    %   Saair Quaderi (off-by-one bug fix, refactoring, parameter changes)
    %

    if nargin < 3
        blurSigmaWidth_pixels = 2;
    end
    hSize = squareSmoothingWindowLen_pixels.*[1 1];
    gaussianKernel = fspecial('gaussian', hSize, blurSigmaWidth_pixels);

    padR = squareSmoothingWindowLen_pixels;
    padC = squareSmoothingWindowLen_pixels;

    blurredImg = padarray(img, [padR, padC], 'replicate', 'both');

    blurredImg = conv2(blurredImg, gaussianKernel, 'same');

    blurredImg = blurredImg((padR + 1):(end - padR), (padC + 1):(end - padC)); 
end
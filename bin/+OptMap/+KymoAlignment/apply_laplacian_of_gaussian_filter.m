function [filteredImg] = apply_laplacian_of_gaussian_filter(img, hSize_pixels, blurSigmaWidth_pixels)
    % APPLY_LAPLACIAN_OF_GAUSSIAN_FILTER - performs a laplacian of Gaussian filter
    %
    % Inputs:
    %   img
    %     the input image
    %   hSize_pixels
    %     the size of the filter's kernel in pixels [height, width]
    %   blurSigmaWidth_pixels
    %     the length of a single standard deviation of the normal
    %     distribution with which the image will be blurred
    %    (in pixels)
    %
    % Outputs:
    %   rtnArr
    %     the filtered image
    %
    % Authors:
    %   Charleston Noble
    %   Saair Quaderi (fixed off-by-one bugs; refactoring)

    padSize = round(hSize_pixels ./ 2);
    imgSize = size(img);

    hSize_pixels = hSize_pixels(1:2);
    logKernel = fspecial('LoG', hSize_pixels, blurSigmaWidth_pixels);
    filteredImg = padarray(img, padSize, 'replicate');

    filteredImg = conv2(filteredImg, logKernel, 'same');
    filteredImg = filteredImg(padSize(1) + (1:imgSize(1)), padSize(2) + (1:imgSize(2)));
end
function [probPostPSF] = apply_point_spread_function(probSeq, psfSigmaWidth, isLinearTF, widthSigmasFromMean)
    % APPLY_POINT_SPREAD_FUNCTION - Performs convolution of a
    %   one-dimensional probability curve with a one-dimensional gaussian
    %   filter kernel representing the camera's point spread function
    %
    % Inputs:
    %   probSeq
    %     Probability values from 0 to 1 for the presence of fluorescence
    %     coming from the one-dimensional position
    %     (e.g. probability of YOYO1 binding to a given basepair along a
    %     linearized DNA sequence)
    %   psfSigmaWidth
    %     The width of a single sigma (mean to a standard deviation from
    %     the mean) of the one-dimensional gaussian distribution that the
    %     probSeq is to beconvolved with.
    %     The unit-length is defined as equivalent to the positional
    %     difference between adjacent values in the probSeq
    %   isLinearTF (optional, defaults to false; e.g. circular convolution)
    %     Whether the convolution should be linear (with values before or
    %     after the sequence being set to 0) as opposed to circular (where  
    %     probability values before or after the sequence are wrapped
    %     around from the opposite end of the sequence)
    %   widthSigmasFromMean (optional, defaults to 4)
    %     the number of deviations from the mean at which the point spread
    %     function's normal distribution is cut to zero in the convolution
    %     kernel
    %
    % Outputs:
    %   probPostPSF
    %     the convolution of the probability sequence with the point spread
    %     function
    %
    % Authors:
    %   Saair Quaderi
    
    if (nargin < 3) || isempty(isLinearTF)
        isLinearTF = false;
    end
    if (nargin < 4) || isempty(widthSigmasFromMean)
        widthSigmasFromMean = 4; %+/- 4 standard deviations from mean
    end

    % psfSigmaWidth must be at same resolution as
    %  indices of probSeq (e.g. basepair resolution)
    
    %%%%% Perform Gaussian convolution on the theory curve.
    nonzeroKernelLen = widthSigmasFromMean*2*psfSigmaWidth;
    nonzeroKernelLen = round((nonzeroKernelLen - 1)/2)*2 + 1; % round to nearest odd integer

    probSeqLen = length(probSeq);
    if isLinearTF
        % Zero-pad sequence probabilities to prevent wrap-around influence
        % from circular convolution
        probSeq = padarray(probSeq, [0, nonzeroKernelLen]);
    elseif nonzeroKernelLen > probSeqLen
        warning(['Probability sequence was shorter than the non-zero PSF kernel length of ', num2str(nonzeroKernelLen)]);
        nonzeroKernelLen = probSeqLen;
    end

    curveLength = length(probSeq);

    % create point spread function kernel
    kernel = fspecial('gaussian', [1, nonzeroKernelLen], psfSigmaWidth);
    psfKernel = zeros(size(probSeq));
    psfKernel(ceil((curveLength - nonzeroKernelLen)/2) + (1:nonzeroKernelLen)) = kernel;

    % shift 50% so that the first and last index are at the peak of the gaussian
    psfKernel = fftshift(psfKernel);

    % compute results from circular convolution
    probPostPSF = cconv(probSeq, psfKernel, curveLength);

    if isLinearTF
        % Extract results out of padded results
        probPostPSF = probPostPSF(nonzeroKernelLen + (1:probSeqLen));
    end
end
function [thyCurve_pxRes] = convert_bpRes_to_pxRes(thyCurve_bpRes, meanBpExt_pixels, islinear)
    % Convert from bp resolution to pixel resolution
    % with moving average window.
    
    if nargin < 3
        islinear = 0;
    end

    if nargin < 2
        s = 541.28;
    else
        s = 1/(meanBpExt_pixels);
    end
    
    
    thyCurve_pxRes(round(length(thyCurve_bpRes)/s)) = 0;
    if islinear
        xtraseq = [zeros(round(s),1); thyCurve_bpRes; zeros(round(s),1)];
    else
        xtraseq = cat(find(size(thyCurve_bpRes) - 1), ...
                      thyCurve_bpRes(end-round(s):end), ...
                      thyCurve_bpRes, ...
                      thyCurve_bpRes(1:round(2*s)));
    end
    for i = 1:round(length(thyCurve_bpRes)/s)
%         data2= xtraseq(floor((i*s)-s+1):floor((i*s)+s));
        thyCurve_pxRes(i) = mean(xtraseq(floor((i*s)-s+1):floor((i*s)+s)));
    end
% 
%     thyCurve_pxRes =[];
%     thyCurve_pxRes(floor(length(thyCurve_bpRes)/s)) = 0;
%     if islinear
%         xtraseq = [zeros(1,round(s)) thyCurve_bpRes zeros(1,round(s))];
%     else
%         xtraseq = [thyCurve_bpRes(end-round(s)+1:end) thyCurve_bpRes thyCurve_bpRes(1:round(s))];
%     end
%     
%     for i = 1:floor(length(thyCurve_bpRes)/(s))
%         thyCurve_pxRes(i) = mean(xtraseq(floor(((i)*s)-s+1):floor(((i*s)))));
%     end

end

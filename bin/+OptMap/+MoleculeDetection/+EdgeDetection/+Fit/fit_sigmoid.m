function [mdl, gofStruct, outputStruct,Fopts] = fit_sigmoid(data,ft,Fopts)
    % fit_sigmoid

    if nargin < 2
        ft = fittype('a + b*normcdf(x,mu,sig)','indep','x');
    end
    
    if nargin < 3
        % fitting options
        Fopts = fitoptions(ft);
        % a, b, mu, sig
        % a - basically want to find mean intensity of background. Lower bound is
        % min(barcode), upper is mean(barcode)
        % b - hight of sigmoid function - distance between the two means. 
        % meanof sigmoid
        % signa of sigmoid - this related to psf 
        % Fopts.Algorithm = 'Levenberg-Marquardt';
        Fopts.Lower = [min(data) 1 0 1];
        Fopts.Upper = [mean(data) max(data)-min(data) length(data) 6];
        Fopts.StartPoint = [1  (max(data)-min(data))/2 10 2.3];

    end

 

	dataPoints = (1:length(data))';

    [mdl, gofStruct, outputStruct]  = fit(dataPoints,data',ft,Fopts);
    


end


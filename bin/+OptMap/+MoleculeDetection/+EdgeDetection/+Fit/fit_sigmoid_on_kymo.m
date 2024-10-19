function [fitsLeft, fitsRight,confLeft,confRight] = fit_sigmoid_on_kymo(kymo)
    % fit_sigmoid_on_kymo 
    % Follows https://se.mathworks.com/matlabcentral/answers/467203-how-to-fit-a-curve-to-a-step-function

    %     Args:
    %        data: kymo
    % 
    %     Returns:
    %       fitsLeft, fitsRight
    import OptMap.MoleculeDetection.EdgeDetection.Fit.fit_sigmoid;

    % initialize
    ft = fittype('a + b*normcdf(x,mu,sig)','indep','x');


    fitsLeft = zeros(size(kymo,1),4);
    fitsRight = zeros(size(kymo,1),4);
    confLeft = zeros(size(kymo,1),2);
    confRight = zeros(size(kymo,1),2);

    % first row
    dataFirst = double(kymo(1,:));
    % deal with nan's, assign first nonnan value to all
    [a,b] = find(~isnan(dataFirst),1);
    dataFirst(isnan(dataFirst)) = dataFirst(b);


    [mdlF, gofStruct, outputStruct, Fopts] = fit_sigmoid(dataFirst,ft);
    fitsLeft(1,:) = [mdlF.a, mdlF.b,mdlF.mu,mdlF.sig];
    % confint
    tempConf = confint(mdlF);
    confLeft(1,:) = tempConf(:,3);
    
    % right side
    Fopts.StartPoint = fitsLeft(1,:) ;
    [mdlR, gofStruct, outputStruct] = fit_sigmoid(fliplr(dataFirst),ft,Fopts);
    fitsRight(1,:) = [mdlR.a, mdlR.b,mdlR.mu,mdlR.sig];
    % conf int
    tempConf = confint(mdlR);
    confRight(1,:) = tempConf(:,3);
    
    
    for k=2:size(kymo,1)
        
        % current row
        dataRow = double(kymo(k,:));
        [a,b] = find(~isnan(dataRow),1);
        dataRow(isnan(dataRow)) = dataRow(b);
        
        % left
        Fopts.StartPoint =  fitsLeft(k-1,:);
        [mdlF, ~, ~, ~] = fit_sigmoid(dataRow,ft,Fopts);
        fitsLeft(k,:) = [mdlF.a, mdlF.b,mdlF.mu,mdlF.sig];
        %conf int
        tempConf = confint(mdlF);
        confLeft(k,:) = tempConf(:,3);
        % right
        Fopts.StartPoint =  fitsRight(k-1,:);
        [mdlR, ~, ~, ~] = fit_sigmoid(fliplr(dataRow),ft,Fopts);

        fitsRight(k,:) = [mdlR.a, mdlR.b,mdlR.mu,mdlR.sig];
        tempConf = confint(mdlR);
        confRight(k,:) = tempConf(:,3);
        
        
    end
    % since this fits the wrong way, value 1 corresponds to size(kymo,2).
    fitsRight(:,3) = size(kymo,2)-fitsRight(:,3)+1;
    confRight(:) = size(kymo,2) - confRight(:)+1;
end


function [ kernel ] = gaussian_kernel( n, sigma )
% 28/11/16
    %kernel = zeros(1,n);
    
    if mod(n,2) == 0
        k = -floor(n/2)+1:1:floor(n/2);
    else
         k = -floor(n/2):1:floor(n/2);
    end
    
    kernel = fftshift(exp(-k.^2./(2*sigma.^2)));
    
    kernel = kernel/sum(kernel);

end


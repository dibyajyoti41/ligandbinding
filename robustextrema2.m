function [ xMinima , xMaxima ] = robustextrema2(curve, IntenThreshold)

    % ROBUSTEXTREMA - Find robust maxima and minuma in a curve
    %
    % Inputs: 
    % curve = 1D curve to be analyzed
    % noiseThreshold = threshold for deeming a peak/valley "robust" or not.
    %
    % Outputs: 
    % xMinima = positions for the robust minima
    % xMaxima = positions for the robust maxima 
    %
    % Dependencies (matlab-functions, MAT-files): none
    %
    % By: Tobias Ambjornsson
    %
    % Refs: 
    % Noble, Charleston, Adam N. Nilsson, 
    % Camilla Freitag, Jason P. Beech, Jonas O. Tegenfeldt, 
    % and Tobias Ambjrnsson. 
    % "A fast and scalable kymograph alignment algorithm for 
    % nanochannel-based optical DNA mappings." 
    % PloS one 10, no. 4 (2015): e0121905, see Supplementary.

  
    diffCurve = diff(curve);
    max_no_boundaries = length(curve);

    x=0;                            % x labels different positions
    boundary_counter=1;             % Any time series can be separated into 
                                    % downward-going and 
                                    % upwards-going regions. 
                                    % The boundary between two such regions 
                                    % is a local minina or local maxima. 
                                    % originates
                                    % from the the DNA melting community
                                    % [a series of papers by Azbel].
    du=zeros(max_no_boundaries,1);         
                                    % positions of down-up left boundaries
 
    ud=zeros(max_no_boundaries,1);         
                                    % positions of up-down left boundaries
    no_of_du_boundaries=0;          % number of du boundaries
    no_of_ud_boundaries=0;          % number of ud boundaries

    M = length(diffCurve);
    while x<=M-1 

       % -- helical segment --

       j=x+1;  % (a possible) energy minimum at position j % first pos is j=1
       C=0;
       while x<=M-1 && C<=IntenThreshold
          x=x+1;
          C=C+diffCurve(x);
          if C<0  % a new possible robust local minima found
             j=x;
             C=0;    
          else % E>=0
             if C>IntenThreshold
                du(boundary_counter)=j;       % store position of the
                                              % start of going up
                no_of_du_boundaries=no_of_du_boundaries+1;  

             end
          end
       end

       % -- coil segment --

       j=x;   % (a possible) energy maximum at position j
       C=0;
       while x<=M-1 && C>=-IntenThreshold
          x=x+1;
          C=C+diffCurve(x);
          if C>0 % a new possible robust local maxima found
             j=x;
             C=0;
          else % E<0
             if C<-IntenThreshold
                ud(boundary_counter)=j;       % store position of
                                              % going down
                no_of_ud_boundaries=no_of_ud_boundaries+1;  

             end
          end
       end

       boundary_counter=boundary_counter+1;    
    end
    
    xMinima = du(du>0)+1;
    xMaxima = ud(ud>0)+1;


    end
    

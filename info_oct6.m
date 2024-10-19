%This code presents the non-GUI version of caculating self information 
%score from theoretical intensity profile or barcode generated from competitive binding
% of two ligands YOYO-1 and Netropsin on \lambda DNA. 
% At first, the ligand arrangements (YOYO-1) are utilized to generate fluoroscence image of DNA.
% The parameters for generating synthetic image are extracted from
% EMCCD-PIA paper by Jens Krog, et. Al, 

% Krog J, Dvirnas A, Strm OE, Beech JP, Tegenfeldt JO, Mu ller V, et al. (2024) 
% Photophysical image analysis: Unsupervised probabilistic thresholding for images from electron-multiplying charge-coupled devices. 
% PLoS ONE 19(4): e0300122. https://doi.org/10.1371/journal. pone.0300122

% Once, multiple time frame images are formed, then kymograph generated ,
% and from kymographs after time-averaging barcode extracted. 

% self information score is calculated on the barcode (intensity profile)
% for a given Netropsin concentration. Self information score checks
% valleys and peaks in the intensity profile in comparison to background
% noise variance introduced originally by
% Noble C, Nilsson AN, Freitag C, Beech JP, Tegenfeldt JO, Ambjrnsson T (2015) 
% A Fast and Scalable Kymograph Alignment Algorithm for Nanochannel-Based Optical DNA Mappings. 
% PLoS ONE 10(4): e0121905. doi:10.1371/journal. pone.0121905
%see equation 4 in the paper.



tic
       noc=0.0001; % netropsin concentration initial choice
      for ji=1:10 % outer loop where netropsin concentration is varied in an interval
       
      for jk=1:500 %sample size 
          % inner loop to generate average self-information score 
          %change it according to the problem. we have used 500 random
          %averaging of information score in the text (sample size 500)
          % (1) from random ligand map \rightarrow (2) image formation (EMCCD PIA) 
          % \rightarrow (3) kymograph and then barcode extracted (4)
          % calculate self-information(jk) score for a given netropsin
          % concentration


       
        numLigands = 2; %competitive binding of 2 ligands
        bindingConstants = [];
        ligandLengths=[4 4]; %Ligand lengths are of four lattice sizes

     
         conc=[noc 0.02];
             
        sets.theoryGen.concDNA=0.2;

        lambdaSequence = fastaread(strcat('sequence.fasta')); % fasta sequence

        tic
        disp('Computing free concentrations');
        ntIntSeq = nt2int( lambdaSequence.Sequence, 'ACGTOnly',1);


      nSize = size(ntIntSeq,2);
      ntNetrospinSeq = zeros(numLigands,nSize);


    
           for i = 1:1
             
             formatSpec = ['%', num2str(ligandLengths(i)), 'c %f'];
             file='binding_constant_rules.txt';
             fid = fopen(file); 
            value = textscan(fid,formatSpec,'delimiter','\n'); fclose(fid);     
            
            ConstantNames = value{1};
            
            ConstantVals = value{2};

            bindingConstantsVals{i}=ConstantVals;
            bindingConstantsNames{i}=ConstantNames;


            seqSpecLen = ligandLengths(i);

           
            bindingConstantsMatSize = repmat(4, [1, seqSpecLen]);
            
            bindingConstantsMat = ones(bindingConstantsMatSize);
            numRules = size(bindingConstantsVals{i}, 1);



            % convert vector to int
            bitsmartTranslationArr = uint8(pow2(3:-1:0));

            for ruleNum=1:numRules
                vect_uint8 = bitsmartTranslationArr(nt2int(bindingConstantsNames{i}(ruleNum,:)));
                mat_logical = logical(rem(floor(double(vect_uint8(:))*pow2(-3:0)),2));
                idxs = mat2cell(mat_logical, ones([1, size(mat_logical, 1)]), 4);
                bindingConstantsMat(idxs{:}) = bindingConstantsVals{i}(ruleNum);
            end

            
            netropsinBindingConstant = 0.4*bindingConstantsMat./1E6;
            


       numbers=0:ligandLengths(i)-1;  
     for j=nSize-ligandLengths(i):-1:1 %1 = 'A' %2 = 'C', 3 = 'G', 4 = 'T'
        str = sprintf('ntIntSeq(%d+j),', numbers); 
        str = str(1:end-1);
       
      
        ntNetrospinSeq(i,j)=eval(['netropsinBindingConstant(', str, ')']);
     end

           end
           
           bindingConstant=[26]; %Binding constant of YOYO-1
         
    for j=2:numLigands
 
          ConstantVals = bindingConstant(j-1);

         bindingConstantsVals{j} = ConstantVals;
        
         ntNetrospinSeq(j,:)=ConstantVals;
    end
    
    % ...


    
        sets.theoryGen.concN=conc(1);
        sets.theoryGen.concY=conc(2);

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=1;
          numbers=0:ligandLengths(i)-1;  
     for j=nSize-ligandLengths(i):-1:1 %1 = 'A' %2 = 'C', 3 = 'G', 4 = 'T'
        str = sprintf('ntIntSeq(%d+j),', numbers); 
        str = str(1:end-1);

        ntNetrospinSeq(i,j)=eval(['netropsinBindingConstant(', str, ')'])*conc(i); %statistical weight of netropsin bound c*k(i)
      end
          
          
    for j=2:numLigands


          ConstantVals = bindingConstant(j-1);

         bindingConstantsVals{j} = ConstantVals;
         ntNetrospinSeq(j,:)=ConstantVals*conc(j); %statistical weight of Yoyo bound c(yoyo)*k_{yoyo}

    end
    
    % ...


  
%%%%%%%%%%%%%%%%%Co-operativity-matrix%%%%%%%%%%%%%%%%%%
% The co-operativity matrix with size (numLigands,numLigands)
% All the entries of the matrix parameters are here set to 1.
 n1=numLigands;
coopc=zeros(1,n1);
coop = zeros(n1);
for i = 1:n1
    for j = 1:n1
        coop(i, j) = 1;
        coopc(i)=1;
    end
end






 len=nSize;



nrun=1; %Simulation run to make average
res=zeros(nrun,len);

s=numLigands;
se=ones(1,s);
sd=ligandLengths;
s1=sum(ligandLengths)+1;
avgprob=zeros(s+1,len);
oyo=zeros(s+1,nrun,len);


marker=zeros(1,s1); 
marker(1)=1;
i1=1;
for i=1:s
    marker(i1+1:i1+sd(i))=sd(i);
i1=i1+sd(i);
end



tic
for k=1:nrun

yoyoBinding = zeros(nSize,1);
leftVec = zeros(nSize+1,s1);     % UL vectors
rightVec = zeros(s1,nSize+1);    %UR vectors
maxEltLeft = zeros(1, nSize+1);
maxEltRight = zeros(1, nSize+1);
U = zeros(s1,s1,s1);
for i=1:s1
    U(i,i,i)=1;
end


leftVecNext(1,:)=zeros(1,s1);
rightVecNext(:,1)=zeros(s1,1);


rightVec(:,nSize+1) = transpose([1 zeros(1,s1-1)]);
maxEltRight(nSize+1) = norm(rightVec(:,nSize+1));

rightVec(:,nSize+1)  = rightVec(:,nSize+1) ./norm( rightVec(:,nSize+1) );   %v(N+1)/|v(N+1)|
ud=ones(1,s1); % the z state vector which determines in which state the site is

    maxVecDiv =  zeros(1,nSize);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%v0%%%%%%%%%%%%%%%%%%%%%%%%
leftVecZero=zeros(1,s1);
leftVecZero(1,1)=1;
leftVecNext=zeros(1,s1);
leftVecPrev = leftVecZero;
leftVecNext(1)=leftVecPrev(1)+leftVecPrev(2)+sum(leftVecPrev(2+cumsum(sd(1:end-1))));
i=0;


b1=0;
b11=1;
for jj=1:s
leftVecNext(sd(jj)+b1)=leftVecPrev(sd(jj)+b1+1)*0;
%%%%% 
b111=0;
 for jjj=1:s
 leftVecNext(sd(jj)+b11)=leftVecNext(sd(jj)+b11)+leftVecPrev(2+b111)*coop(jjj,jj);
 b111=b111+sd(jjj);
 end
  leftVecNext(sd(jj)+b11)=leftVecNext(sd(jj)+b11)+leftVecPrev(1);
  b11=b11+sd(jj);

b1=b1+sd(jj);
end

j=1;
j1=1;
while j<s+1
      bb1=cumsum(sd(1:j));
      for jj=j1+1:bb1(end)-1
     leftVecNext(jj)=leftVecPrev(jj+1);
      end
      j1=bb1(end)+1;
j=j+1;
end 

maxEltLeftZero = norm(leftVecNext); %nL
leftVec(i+1,:) = leftVecNext./maxEltLeftZero; %wL=uL/nL

%%%%%%%%%%%%%%%%%%%%%%%%%%v0%%%%%%%%%%%%%%%%%%%%%%%%




%====================================================
i=1;
stateChosen=ones(1,len); %the states m_i
PCond= zeros(len,s1);


for i=1:nSize
leftVecNext=zeros(1,s1);
leftVecPrev = leftVec(i,:);
leftVecNext(1)=leftVecPrev(1)+leftVecPrev(2)+sum(leftVecPrev(2+cumsum(sd(1:end-1))));




b1=0;
b11=1;
for jj=1:s
leftVecNext(sd(jj)+b1)=leftVecPrev(sd(jj)+b1+1)*ntNetrospinSeq(jj,i);
%%%%% 
b111=0;
 for jjj=1:s
 leftVecNext(sd(jj)+b11)=leftVecNext(sd(jj)+b11)+leftVecPrev(2+b111)*coop(jjj,jj);
 b111=b111+sd(jjj);
 end
  leftVecNext(sd(jj)+b11)=leftVecNext(sd(jj)+b11)+leftVecPrev(1);
  b11=b11+sd(jj);

b1=b1+sd(jj);
end

j=1;
j1=1;
while j<s+1
      bb1=cumsum(sd(1:j));
      for jj=j1+1:bb1(end)-1
     leftVecNext(jj)=leftVecPrev(jj+1);
      end
      j1=bb1(end)+1;
j=j+1;
end 

maxEltLeft(i) = norm(leftVecNext); %nL
leftVec(i+1,:) = leftVecNext./maxEltLeft(i); %wL=uL/nL
end        

    ud=eye(s1);
    stateSeq = zeros(1,nSize);
   
    changeStates= [1 2 2+cumsum(sd(1:end-1))];
    relVectors = ud(changeStates,:);
    nextStates = [];
    i=1;
% while i<nSize+1
 for i=1:nSize
         rightVecPrev = rightVec(:,nSize-i+2);
         b2=1;
        
       % rightVecPrev=rightVecPrev2'*ud;

         rightVecNext=zeros(s1,1);
         rightVecNext(1)=rightVecPrev(1);
         b21=2;
         for jj=1:s
               % b21=2;
            rightVecNext(1)=0;
             b22=1;
             for jj1=1:s
                 b22=b22+sd(jj1);
                 rightVecNext(1)=rightVecNext(1)+rightVecPrev(b22);
                 rightVecNext(b21)=rightVecNext(b21)+rightVecPrev(b22)*coop(jj,jj1);       
             end
             rightVecNext(b21)=rightVecNext(b21)+rightVecPrev(1);
             rightVecNext(1)=rightVecNext(1)+rightVecPrev(1);
             
             b21=b21+sd(jj);
         end
     

j1=2;
for j=1:s
      bb1=cumsum(sd(1:j));
      for jj=j1+1:bb1(end)
      rightVecNext(jj)=rightVecPrev(jj-1);
      end
      rightVecNext(bb1(end)+1)=rightVecPrev(bb1(end))*ntNetrospinSeq(j,nSize-i+1);
      j1=bb1(end)+2;

end 
       
          %rightVecNext=ud*rightVecNext;
         if ~isempty(nextStates)
            stateSeq(i) = nextStates(1);
            nextStates(1) = [];
            rightVecNext = ud( stateSeq(i),:)'.*rightVecNext;
            randNum = 1;
         else
            rightVecUpdated = relVectors'.*rightVecNext;
            probs = leftVec(nSize-i+1,:)*rightVecUpdated/(leftVec(nSize-i+2,:)*rightVecPrev)/maxEltLeft(nSize-i+1);
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        a1=probs(probs>0);
        x=rand;  
        j=1;
        k2=length(a1);
        while j<k2

        mid=floor((j+k2)/2);
        p=cumsum(a1(1:mid));
        if x<p(mid)
            k2=mid;
        else
            j=mid+1;
        end
        end 

         e1=find(probs==a1(j));
          stateSeq(i)=changeStates(e1);
          rightVecNext = rightVecUpdated(:,e1);
          randNum=stateSeq(i);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         end
          maxEltRight(nSize+1-i) = norm(rightVecNext);
        rightVec(:,nSize+1-i) = rightVecNext./maxEltRight(nSize+1-i);

        if (randNum > 1) % then need to fast forward to the next states
            nextStates = stateSeq(i)+1:1:stateSeq(i)+(marker(stateSeq(i))-1);
        end
    
 
         
end
stateChosen=stateSeq;
res(k,:)=stateChosen;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
te=([]);
te = find(stateChosen == 1);


for i = 1:s
    if i == 1
        temp1 = find(stateChosen > 1 & stateChosen < 1 + sd(i) + 1);
        te=[te; {temp1}];
    else
     
     temp2=find(stateChosen > sum(sd(1:i - 1)) + 1 & stateChosen < 2 + sum(sd(1:i)));
    te=[te; {temp2}];
    end

end




for i=1:s+1
    oyo(i,k,te{i,1})=1;
end 




%%%%%%%%%%%%%%calculatingyoyosineachpixel%%%%%%%%%%%%
sets.theoryGen.pixelWidth_nm=110; % pixel size 110 nm
sets.theoryGen.meanBpExt_nm=0.3;  %nmBp ratio


numLoops = floor(nSize/(sets.theoryGen.pixelWidth_nm/sets.theoryGen.meanBpExt_nm)); %length of DNA in pixels
entriesPerLoop = floor(sets.theoryGen.pixelWidth_nm/sets.theoryGen.meanBpExt_nm); % number of bps in each pixel
for i = 1:numLoops
startIndex = (i - 1) * entriesPerLoop + 1;
endIndex = startIndex + entriesPerLoop - 1;
Inten(k,i)= sum(oyo(3,k,startIndex:endIndex)); %number of yoyos in each pixel 
end


%%%%%%%%%%%%%%calculatingyoyosineachpixel%%%%%%%%%%%%

end
Inten=Inten/4;

toc


%%%%%%%%%%%%%%%%cameranoiseparameters(sCMOS)%%%%%%%%%%%%
  %% parameters are taken from EMCCD-PIA table 
 
  %Krog J, Dvirnas A, Stro m OE, Beech JP, Tegenfeldt JO, 
  % Mu ller V, et al. (2024) 
  % Photophysical image analysis: Unsupervised probabilistic thresholding for images from electron-multiplying charge-coupled devices. 
  % PLoS ONE 19(4): e0300122. https://doi.org/10.1371/journal. pone.0300122

    chipPars.NA  = 1.49;
    chipPars.wavelength = 660;
    chipPars.pixelsize = 110;
    chipPars.gain = 18.82; % EM-Gain
    chipPars.adfactor = 35.17; % AD factor 1/f
    chipPars.countoffset = 26.37;%offset
    chipPars.ronoise = 1.45;%r0 noise
    chipPars.QE = 1.0; %Quantum efficiency (0.95)
    chipPars.c = 0.002; %constant

    snr=10;
    noisePars.lambdaBg = 22;  % Poisson parameter for background (\lambda_{Bg})
    %noisePars.lambdaSg = 28.4; % Poisson parameter for signal (\lambda_{Sig})
    %noisePars.lambdaSg = 0.5*snr^2*(1+(1+4*noisePars.lambdaBg/snr^2)^0.5) ; % Poisson parameter for signal (\lambda_{Sig})
    noisePars.pxSide = 20; 
    %noisePars.sigma = 1.845;   % Background noise variance 
    noisePars.sigma = 300/110; %\sigma_{psf}  

    
    
    labelPars.nmbp = sets.theoryGen.meanBpExt_nm;      % nanometer-basepair ratio;
    labelPars.frameLen = 0.1;     
    labelPars.numFrames=50; %number of timeframes in synthetic image
    labelPars.fTime = labelPars.frameLen*labelPars.numFrames;
    
%     eyo=numLoops*noisePars.lambdaSg/sum(Inten);

    
    nX = 256; % size of the image
    nY = 32; % fluorophores will be centered at nY/2
  
    ngapPeaks = 15000;
   
        dpeaks = floor(nSize/(numLoops));
       
   
        npeaks = numLoops;
    
        blinking = 0;
    
    
    [fluorophorePos] = fixed_fluorophores(npeaks,ngapPeaks,dpeaks,Inten);
    fluorophorePosPx = [fluorophorePos*labelPars.nmbp/chipPars.pixelsize round(nY/2)*ones(length(fluorophorePos),1)];
    
   
    [fluorophorePosMulti] = label_position_sim(fluorophorePos',labelPars);% positions,totLen,numFrames,frameLen);

    
        photonProb = ones(size(fluorophorePosMulti));
    
    final1D = zeros(size(photonProb,2),nX);
    finalImage = cell(1,size(photonProb,2));
  
   parfor i=1:size(photonProb,2)
        fluorophorePosPxCur = [fluorophorePosMulti(:,i)*10^9/chipPars.pixelsize round(nY/2)*ones(length(fluorophorePos),1)];
        
        Intenn=Inten;
        [pxPhotons, img] = gen_photon_prob(Intenn,noisePars.lambdaBg,fluorophorePosPxCur,noisePars.sigma, nX, nY);
        
        noisyImage=noise_model(pxPhotons',chipPars.gain,chipPars.ronoise,chipPars.adfactor,chipPars.countoffset,chipPars.c,chipPars.QE,noisePars.lambdaBg);
       
        finalImage{i} = reshape(noisyImage,nY,nX); %this is the formed image with dimension nx \times ny =256 \times 32 (fig 6 (a), (d) and (g))

        B=squeeze(finalImage{i});
       
        final1D(i,:) = finalImage{i}(round(nY/2),:);
    end
    %figure,imagesc(finalImage{1}) %uncomment this to generate the image
%%%%%%%%%%%%%%%%cameranoiseparameters%%%%%%%%%%%%

dbmStruct.kymoCells.rawKymos = cell(1,1);
dbmStruct.kymoCells.rawKymos{1} = zeros(labelPars.numFrames,nX);
dbmStruct.kymoCells.rawBitmask{1} = zeros(labelPars.numFrames,nX);
for i=1:labelPars.numFrames

dbmStruct.kymoCells.rawKymos{1}(i,:) = mean(finalImage{1,i}(nY/2-2:nY/2+2,:)); %rawkymos generated from the image
dbmStruct.kymoCells.rawBitmask{1}(i,round(fluorophorePosMulti(1,i)*10^9/chipPars.pixelsize)-0:round(fluorophorePosMulti(end,i)*10^9/chipPars.pixelsize)+0)=1;
end

useGUI = 0;



numK=1;
if ~isfield(dbmStruct.kymoCells,'alignedKymos')
        import OptMap2.KymoAlignment2.NRAlign2.nralign; % aligned the kymograph accordingly the features
        %import OptMap.KymoAlignment.NRAlign.nralign;
            for ix=1:numK
                    fprintf('Aligning kymograph for file molecule #%d of #%d ...\n', ix, numK);
                    [alignedKymo, stretchFactorsMat, shiftAlignedKymo,alignedBitmask] = nralign(dbmStruct.kymoCells.rawKymos{ix},false,dbmStruct.kymoCells.rawBitmask{ix});
                    dbmStruct.kymoCells.alignedKymos{ix} = alignedKymo;  %aligned kymograph %this is plotted as in fig 6 (b),(e), and (h)
                    dbmStruct.kymoCells.alignedBitmask{ix}=alignedBitmask;
                    cg=nanmean(dbmStruct.kymoCells.alignedKymos{ix},1);
                    cg2=nanmean(dbmStruct.kymoCells.alignedBitmask{ix},1);
                    dbmStruct.kymoCells.alignedBitmasks{ix} = ~isnan(alignedKymo);
                    dbmStruct.kymoCells.stretchFactorsMat{ix} = stretchFactorsMat;
                    dbmStruct.kymoCells.shiftAlignedKymo{ix} = shiftAlignedKymo;     
                   
            end 
end
           


           for ix=1:numK
                dbmStruct.kymoCells.barcodes{ix} = nanmean(dbmStruct.kymoCells.alignedKymos{ix},1);
                dbmStruct.kymoCells.barcode{ix} = cg(cg2>0); %barcode extracted
                dbmStruct.kymoCells.barcodesStd{ix} = nanstd(dbmStruct.kymoCells.alignedKymos{ix},1,1);
                dbmStruct.kymoCells.numsKymoFrames(ix) = size(dbmStruct.kymoCells.alignedKymos{ix},1);
                try
                    bitmask = dbmStruct.kymoCells.alignedMask{ix};
                catch
                    bitmask = ~isnan(dbmStruct.kymoCells.alignedKymos{ix});
                end
                dbmStruct.kymoCells.fgStartIdxs{ix} = arrayfun(@(x) find(bitmask(x,:),1,'first'),1:size(bitmask,1));
                dbmStruct.kymoCells.fgEndIdxs{ix} = arrayfun(@(x) find(bitmask(x,:),1,'last'),1:size(bitmask,1));
                dbmStruct.kymoCells.fgStartIdxsMean(ix) = round(mean( dbmStruct.kymoCells.fgStartIdxs{ix}));
                dbmStruct.kymoCells.fgEndIdxsMean(ix) = round(mean(  dbmStruct.kymoCells.fgEndIdxs{ix}));
            end
%oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

%figure,imagesc(dbmStruct.kymoCells.alignedKymos{ix})%uncomment this to see
%the kymograph



vectIn=dbmStruct.kymoCells.barcode{1}; 


%%%%%%%%meanvarianceupdated%%%%%%%%%%%%
cf=2*(chipPars.gain/chipPars.adfactor);%commonfactor

sigtot=(cf*mean(vectIn)-(cf*chipPars.countoffset)+(chipPars.ronoise)^2+(1/12))^0.5;
%%%%%%%%meanvarianceupdated%%%%%%%%%%%%

    toc
   
    %%%%%%%%%%%%%%infoscore calculation%%%%%%%%%%%%%%%%%%%%%%%
   % ba=dbmStruct.kymoCells.rawKymos{1,1}(:,200:end); % background noise variance calculation for the image part without DNA
    ba=finalImage{1,1}(:,200:end);
    ba1=ba(:);
   % estNoiseStd=std(ba1);
    noiseThreshold=std(ba1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    energyLandscape=vectIn;
  
    IntenThreshold=6*sigtot/(50^0.5); %I_th=X\sigma_{th}/(num_of_frames)^0.5


    curve=energyLandscape;
    [ xMinima , xMaxima ] = robustextrema2(curve, IntenThreshold);
    



    infoScore = info_score(noiseThreshold,xMinima,xMaxima);
    
    
    
      fid = fopen('scorenet_octx6_updated.txt', 'a+');
      fprintf(fid, '%2.5f %2.2f\n', conc(1), infoScore);
      fclose(fid);

      toc
      
      clear
       aab=load('scorenet_octx6_updated.txt');
      noc=aab(end,1);

      end
      
      aab=load('scorenet_octx6_updated.txt');

      noc=aab(end,1)*10;
      
      end

%%%%%%%%%%%%%%%%%result part (Fig 5 in the main text)%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%plot average self information score w netropsin concentration%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



ab1=load('scorenet_octx6_updated.txt');

hk=0;

wnd=500;

for i=1:10 %%%this is for the no of concentrations for netropsin
    ab2=ab1(hk+1:hk+wnd,2);
    mh(i)=mean(ab2(~isinf(ab2)));
    
  
    mh3(i)=std(ab2(~isinf(ab2)));
    hk=hk+wnd;
end




x(1)=0.0001;
for i=1:9
    x(i+1)=x(i)*10;
   
end
figure('units','normalized','outerposition',[0 0 1 1])

set(0,'DefaultLineLineWidth',2)




mean_data=mh;

std_dev=mh3;

upper_bound = mean_data + 1 * std_dev;
lower_bound = mean_data - 1 * std_dev;

h = errorbar(x, mean_data, std_dev, '-o'); % Create the plot with default colors

% Set the line and marker color to blue
h.Color = 'b';
h.MarkerEdgeColor = 'b';
h.MarkerFaceColor = 'b';


h.LineWidth = 1.5; 
h.CapSize = 10;


errorbarLineHandles = get(h, 'Children');
% Find the line objects that correspond to the error bars
for k = 1:length(errorbarLineHandles)
    if strcmp(errorbarLineHandles(k).Type, 'line')
        set(errorbarLineHandles(k), 'Color', 'r'); 
    end
end

% Optional: Add labels and title

hold on

xlabel('Netropsin conc')
ylabel('Infoscore')
legend('Score with errorbars', 'Simulation(mean)','Simulation(mean)\pm\sigma')
legend('X=6')












%%%%%%%%%%%%%%%%%%equidistantfluorophoresalongdna%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fluorophorePos] = fixed_fluorophores(npeaks,ngapPeaks,dpeaks,Inten)

firstPos = ngapPeaks:ngapPeaks:npeaks*ngapPeaks;
secondPos=  firstPos+dpeaks; % dpeaks tunable

fluorophorePos = ngapPeaks:dpeaks:ngapPeaks+(dpeaks*(length(Inten)-1));

fluorophorePos=fluorophorePos';
end
%%%%%%%%%%%%%%%%%%GaussianPSF%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [pxPhotons,img] = gen_photon_prob(Inten2,eeyo,lambdas,lambdaSig,lambdaBg,fluorophorePosPx,sigma,nX,nY)
function [pxPhotons,img] = gen_photon_prob(Inten2,lambdaBg,fluorophorePosPx,sigma,nX,nY)

[X,Y] = meshgrid(1:nX,1:nY);

X = X(:);
Y = Y(:);
pxPhotons = zeros(size(X));

deltaEx = @(x,b) 1/2*(erf((x-b+1/2)/(sqrt(2)*sigma))-erf((x-b-1/2)/(sqrt(2)*sigma)));
j1=1;
for i=1:size(fluorophorePosPx,1)
     %pxPhotons = pxPhotons +4*eeyo*Inten2(i)*deltaEx(X,fluorophorePosPx(i,1)).*deltaEx(Y,fluorophorePosPx(i,2));
     pxPhotons = pxPhotons +25*Inten2(i)*deltaEx(X,fluorophorePosPx(i,1)).*deltaEx(Y,fluorophorePosPx(i,2)); %25=\lambda_{sig}/YOYO
    

end

pxPhotons = pxPhotons+lambdaBg;
img = reshape(pxPhotons,nY,nX);


end


%%%%%%%%%%%%%%%%%%%%EMCCD-PIA:synthetic image generation%%%%%%%%%%%%%%%%%
function noisyImage=noise_model(photonsPlaced,gain,roNoise,adFactor,offset,c,QE,lambdaBg)

    poissVals = poissrnd(photonsPlaced*QE+c);
    if gain == 1
        elecImage = poissVals;
    else
        elecImage = gamrnd(poissVals,gain);
    end
   
    noisyImage = (roNoise) * randn(1,length(photonsPlaced)) + offset + elecImage/adFactor;
    
    
    noisyImage = round(noisyImage);

end


%%%%position of fluorophoroes in multiframes%%%%%%%%%%
function [circPos] = label_position_sim(map,labelPars)

    % parameters in labelPars
    nmbp = labelPars.nmbp;
    numFrames = labelPars.numFrames; % sample frames based on this sampling rate 
   
    pos = map*nmbp/10^9;
   
    circPos(:,1) = pos;
    for i=2:numFrames
          circPos(:,i) = pos;
    end


end









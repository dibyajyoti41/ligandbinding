

%%%In this non-Gui version we start with the defined inputs such as 
     %Input lattice template=  λ- DNA fasta sequence(48502 lattice sites )
     % Number of Ligands=3
     % Ligand Length= 4(ligand 1), 4(ligand 2), and 4(ligand 3)
     %Concentration of Ligand 1 =6 mM, Ligand 2= 0.02 mM and Ligand 3= 0.02 mM
     % Binding Constant of lIgand 1= sequence dependent = Extracted from DL
     % Boger fluroscence intensity chart.
     %Binding constant of Ligand 2= 26x 10^6 and Ligand3=10x 10^6
     % Simulation runs 200 and then the average will be compared with exact
     % transfer matrix calculation which appended in this same code at
     % last.


      
       numLigands = 3; 
        bindingConstants = [];
        ligandLengths=[4 4 4]; %Ligand lengths are of four lattice sizes
        conc=[6 0.02 0.02];  

       
        sets.theoryGen.concDNA=0.2;

        lambdaSequence = fastaread(strcat('sequence.fasta'));

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
       
       % ntNetrospinSeq(i,j)=eval(['netropsinBindingConstant(', str, ')'])*conc(i);
        ntNetrospinSeq(i,j)=eval(['netropsinBindingConstant(', str, ')']);
     end

           end
           
           bindingConstant=[26 10]; %Binding constants 
         
    for j=2:numLigands
 
          ConstantVals = bindingConstant(j-1);

         bindingConstantsVals{j} = ConstantVals;
        
         ntNetrospinSeq(j,:)=ConstantVals;
    end
    
    % ...


    
        sets.theoryGen.concN=conc(1);
        sets.theoryGen.concY=conc(2);
%        sets.theoryGen.concEtbr=conc(3);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=1;
          numbers=0:ligandLengths(i)-1;  
     for j=nSize-ligandLengths(i):-1:1 %1 = 'A' %2 = 'C', 3 = 'G', 4 = 'T'
        str = sprintf('ntIntSeq(%d+j),', numbers); 
        str = str(1:end-1);

        ntNetrospinSeq(i,j)=eval(['netropsinBindingConstant(', str, ')'])*conc(i);
      end
          
          
    for j=2:numLigands


          ConstantVals = bindingConstant(j-1);

         bindingConstantsVals{j} = ConstantVals;
         ntNetrospinSeq(j,:)=ConstantVals*conc(j);

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





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nSize=50;
 len=nSize;



nrun=100; %Simulation run to make average
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

% 
%leftVec(1,:)=zeros(1,s1);
leftVecNext(1,:)=zeros(1,s1);
rightVecNext(:,1)=zeros(s1,1);
%leftVec(1,1+cumsum(sd))=1;
%leftVec(1,1)=1;


%leftVec(1,:) = leftVec(1,:)./norm(leftVec(1,:)); %v(1)/|v(1)|

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
%leftVecNext(1+cumsum(sd))=leftVecNext(1);
% 



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

%leftVecNext(1+cumsum(sd))=leftVecNext(1);
% 



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


end


toc


%%%%%%%%%%%%%%%%%%%%%%%Theorypart%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
leftVect = zeros(nSize+1,s1);
rightVect = zeros(s1,nSize+1);
maxEltLeftt = zeros(1, nSize);
maxEltRightt = zeros(1, nSize);

% 
leftVect(1,:)=zeros(1,s1);
leftVecNextt(1,:)=zeros(1,s1);
rightVecNextt(:,1)=zeros(s1,1);
leftVect(1,1+cumsum(sd))=1;
leftVect(1,1)=1;


leftVect(1,:) = leftVect(1,:)./norm(leftVect(1,:));

rightVect(:,nSize+1) = transpose([1 zeros(1,s1-1)]);


rightVect(:,nSize+1)  = rightVect(:,nSize+1) ./norm( rightVect(:,nSize+1) );
ud=ones(1,s1);

    maxVecDivt =  zeros(1,nSize);


for i=1:nSize
leftVecNextt=zeros(1,s1);
leftVecPrevt = leftVect(i,:);
leftVecNextt(1)=leftVecPrevt(1)+leftVecPrevt(2)+sum(leftVecPrevt(2+cumsum(sd(1:end-1))));
%leftVecNextt(1+cumsum(sd))=leftVecNextt(1);


b1=0;
b11=1;
for jj=1:s
leftVecNextt(sd(jj)+b1)=leftVecPrevt(sd(jj)+b1+1)*ntNetrospinSeq(jj,i);
b111=0;
for jjj=1:s
leftVecNextt(sd(jj)+b11)=leftVecNextt(sd(jj)+b11)+leftVecPrevt(2+b111)*coop(jjj,jj);
b111=b111+sd(jjj);
end
 leftVecNextt(sd(jj)+b11)=leftVecNextt(sd(jj)+b11)+leftVecPrevt(1);
 b11=b11+sd(jj);
 %above
b1=b1+sd(jj);
end

j=1;
j1=1;
while j<s+1
      bb1=cumsum(sd(1:j));
      for jj=j1+1:bb1(end)-1
      leftVecNextt(jj)=leftVecPrevt(jj+1);
      end
      j1=bb1(end)+1;
j=j+1;
end 
maxEltLeftt(i) = norm(leftVecNextt);
leftVect(i+1,:) = leftVecNextt./maxEltLeftt(i);

end



for i=1:nSize
         rightVecPrevt = rightVect(:,nSize-i+2);
         b2=1;

         rightVecNextt=zeros(s1,1);
         rightVecNextt(1)=rightVecPrevt(1);
         b21=2;
         for jj=1:s

            rightVecNextt(1)=0;
             b22=1;
             for jj1=1:s
                 b22=b22+sd(jj1);
                 rightVecNextt(1)=rightVecNextt(1)+rightVecPrevt(b22);
                 rightVecNextt(b21)=rightVecNextt(b21)+rightVecPrevt(b22)*coop(jj,jj1); 
                 
             end
             rightVecNextt(b21)=rightVecNextt(b21)+rightVecPrevt(1);
             rightVecNextt(1)=rightVecNextt(1)+rightVecPrevt(1);

             b21=b21+sd(jj);
         end
j1=2;
for j=1:s
      bb1=cumsum(sd(1:j));
      for jj=j1+1:bb1(end)
      rightVecNextt(jj)=rightVecPrevt(jj-1)*ud(jj-1);
      end
      rightVecNextt(bb1(end)+1)=rightVecPrevt(bb1(end))*ntNetrospinSeq(j,nSize-i+1);
      j1=bb1(end)+2;

end 




         maxEltRightt(nSize+1-i) = norm(rightVecNextt);

         rightVect(:,nSize+1-i) = rightVecNextt./maxEltRightt(nSize+1-i);

end  

    maxVecDivt =  zeros(1,nSize);
    maxVecDivt(1) = maxEltLeftt(1)/maxEltRightt(1);

    for i=2:nSize
        maxVecDivt(i) = maxVecDivt(i-1)*maxEltLeftt(i)/maxEltRightt(i);
    end


    denominator = leftVect(1,:)*rightVect(:,1);




for i=1:s+1
oMat(:,:,i) = zeros(s1, s1);
end
oMat(1,1,1)=1;

cs=1;
for j=2:s+1
for i=1:sd(j-1)
    oMat(i+cs,i+cs,j)=1;
end
cs=cs+sd(j-1);
end

ptheory=zeros(s+1,nSize);
for i=1:s+1
    ptheory(i,1)=leftVect(1,:)*oMat(:,:,i)*rightVect(:,1)/denominator;
 end

for i=1:s+1
    for j=2:nSize

    ptheory(i,j)=leftVect(j,:)*oMat(:,:,i)*rightVect(:,j)*maxVecDivt(j-1)/denominator;
    end
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%






% figure,plot(outputArray4)
%%
% 
%  $PREFORMATTED
%  TEXT$ 
% 
x=1:300;
figure('units','normalized','outerposition',[0 0 1 1])
%figure
set(0,'DefaultLineLineWidth',2)
% 
% First subplot
for i=1:s+1

subplot(s+1,1,i)

m1d = mean(oyo(i,:,:),2);
tempd=squeeze(m1d);
m2d=(tempd)';
mean_data=fliplr(m2d);
s1d = std(oyo(i,:,:));
s2d=squeeze(s1d);
s3d=(s2d)';
std_dev=fliplr(s3d);

upper_bound = mean_data + 1 * std_dev;
lower_bound = mean_data - 1 * std_dev;

plot(1:300,ptheory(i,1:300), 'b-',1:300,mean_data(1:300) , 'r--')

hold on
fill([x, fliplr(x)], [upper_bound(1:300), fliplr(lower_bound(1:300))], 'g', 'EdgeColor', 'k');
alpha(0.3); % Adjust transparency if needed
hold on
if (i==1)
title('Free Sites')
xlabel('Length in Basepairs')
ylabel('Probability')
else
    %sprintf('%dth col of X vs %dth col of Y',k,k)
    title(sprintf('Ligand%d',i))
xlabel('Length in Basepairs')
ylabel('Probability')
legend('Theory', 'Simulation(mean)','Simulation(mean)\pm\sigma')
end

end 





%%%%%%%%%%%%%%%%%%%%%%withsamplingmeaan%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('units','normalized','outerposition',[0 0 1 1])
%figure
set(0,'DefaultLineLineWidth',2)
for i=1:s+1

subplot(s+1,1,i)

ayo2=fliplr(squeeze(oyo(i,:,:)));
ayo=ayo2(:,1:300);
N = 10;% sub-sample to get std
mn2 = arrayfun(@(x) mean(ayo(x:x+N-1,:)),1:N:nrun-N+1,'UniformOutput',false);

 yo1= mean(ayo);
 yo1std = std(cell2mat(mn2')); 



%%
%figure

errorbar(yo1,yo1std,'LineWidth',1.5)
hold on
%plot(flipud(ptheory(3,1:300))) 
plot(1:300,ptheory(i,1:300)) 
hold on
if (i==1)
title('Free Sites')
xlabel('Length in Basepairs')
ylabel('Probability')
legend('Simulation','Theory')
else
    %sprintf('%dth col of X vs %dth col of Y',k,k)
    title(sprintf('Ligand%d',i))
xlabel('Length in Basepairs')
ylabel('Probability')
legend('Simulation','Theory')
end


end

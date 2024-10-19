function [maxcoef, pos, or, secondPos, lenM,dist] = ucr_dtw_score(query, data, querybit, databit, R)
    %   ucr_dtw_score - computes dtw score based on "trillion" code from UCR
    %
    %
    %   Args:
    %       shortVec, longVec, shortVecBit
    %   returns:
    %       rezMax - which stores maxcoef,pos, and or
    %
    
%     if nargin < 4   
%         % Sakoe-Chiba band, this corresponds to stretch factor
% %         R = sets.theory.stretchFactors(end)-1;
%         R = 0.01;
%         nameFiles = num2str(randi(100));
%         dtwscriptpath = fullfile(pwd,'ucr_dtw.sh');
%         matDirpath = fullfile(pwd,'output');
%     end

    % for query, only select the non-nan bit
    query = query(logical(querybit));
    locStart = find(query,1,'first');
    
%     R= 5;  % describes the Sakoe-Chiba band, so in this case we don't allow any local stretching
    % and instead increase this if we want better output
%   
% mex 'NAN_DTW_MEX.cpp';
    % check: what happens if bitmask is all ones?
%     tic
    [posD(1), score1] = NAN_DTW_MEX(data,query,double(databit<50), length(data), length(query), R);
    [posD(2), score2] = NAN_DTW_MEX(data,fliplr(query),double(databit<50), length(data), length(query), R);
%     toc
    
    % take the negative of minscore - then this will be max
    [mincoef,or] = min([score1,score2]);
    maxcoef = -mincoef;
    pos = posD(or)+1; % if we want location for the full query
    
    secondPos = locStart;
    lenM = length(query);

    %     
    % relation of R to MAXSAMP?
    %     seems ok!!
    b1 = zscore(data(pos:pos+length(query)-1),1);
    b2 = zscore(query,1);
    if or==2
        b2=fliplr(b2);
    end

    [~,IX,IY] = dtw(b1,b2,'squared',R);
    
    dist=[IX IY];
    
    
    % DIST
%     sqrt(DIST) % gives the same output!
% 
% %plotting stuff:
% figure
% hold on
% v1=zscore(b1,1)+3;
% plot(v1)
% v2=zscore(b2,1)-3;
% plot(v2)
% for i=1:10:length(IX)
%     plot([IX(i),IY(i)],[v1(IX(i)),v2(IY(i))]);
% end
% figure
% plot(IX,IY,'o-',[IX(1) IX(end)],[IY(1) IY(end)])


    % todo: recover the score and position on the 2 barcodes based on the
    % position found here:
%     figure,plot(zscore(data(pos:pos+length(query)-1),1))
%     hold on
%     plot(zscore(query,1))
    %     
%     [dist,ix,iy] = dtw(zscore(query,1),zscore(data(pos:pos+length(query)-1),1),45);
%     sqrt(dist)

%     
%     [dist,ix,iy] = dtw(zscore(query,1),zscore(data(pos:pos+length(query)-1),1),45);
%     sqrt(dist)


%     [pos(1,1), scores(1,1)] = OVERLAPPING_DTW_MEX(theory,shortVec, length(shortVec), 1.05);
%     [pos(1,2), scores(1,2)] = OVERLAPPING_DTW_MEX(theory,flipud(shortVec), length(shortVec), 1.05);

% 
%     % needs to compute dtw
%     shortVecCut = shortVec(logical(shortVecBit));
% 
%     % rand number, later change this to idx of barcode 
% %     nameFiles = sets.idx;
%     % length of experiment
%     M = sum(shortVecBit);
%     % save experiment in temporary txt file. Check how this behaves in case
    % parfor is used
    
    % should regulate the precision via settings file..
%     
%     fname1 = strcat([nameFiles 'query.txt']); fileID = fopen(fname1,'w');
%     fprintf(fileID,'%2.5f ',shortVecCut); fclose(fileID);
%     
%     
%     fname2 = strcat([nameFiles 'queryrev.txt']); fileID = fopen(fname2,'w');
%     fprintf(fileID,'%2.5f ',fliplr(shortVecCut)); fclose(fileID);
%     
%     % this is a bit dumb since we save the theory again
%     fname3 = strcat([nameFiles 'theory.txt']); fileID = fopen(fname3,'w');
%     fprintf(fileID,'%2.5f ',theory); fclose(fileID);
%     

    
%     pathToScript = fullfile(sets.dtwscriptpath,'ucr_dtw.sh');
% 
%     pathToScript = dtwscriptpath;
% 
%     outFile = fullfile(matDirpath,strcat([nameFiles 'output.txt']));
% %     tic
%     ucrCode = fullfile(strrep(dtwscriptpath,'ucr_dtw.sh','a.out'));
% 
%     cmdStr       = [pathToScript ' ' fname1 ' ' fname2 ' ' fname3 ' ' num2str(M) ' ' num2str(R) ' ' matDirpath ' ' outFile ' ' ucrCode];
%     system(cmdStr);
% %     toc
%     delete(fname1);
%     delete(fname2);
%     delete(fname3);

% 
%     % smart would be to use mpi to save to different parts of the file..
%     A = importdata(outFile);
%     delete(outFile);
%     
%     % just make the coeff negative, so we look for max instead of min
%     coef = -[A(2) A(4)];
%     pos = [A(1) A(3)];
%     [rezMax.maxcoef,rezMax.or ] = max(coef);
%     rezMax.pos = pos(rezMax.or)+1-find(shortVecBit,1,'first')+1;
    
    
%     % mex the c function. This Should be mexed before    
%     mex 'bin/OVERLAPPING_DTW_MEX.cpp';
    
	% Instead of saving these every time, just have them saved as txt's and
	% pass them to ucr function. Should be a more clever way to do this
	% though.

%     % reversed
%     fname = strcat(['experimentrev.txt']); fileID = fopen(fname,'w');
%     fprintf(fileID,'%2.5f ',fliplr(shortVec(shortVecBit))); fclose(fileID);
%     
%     
%     % subsequence length.
%     N = sum(shortVecBit);
% 
%         % these first two are just to see how the query places as a
%     % subsequence, they're not comparable to the rest of the scores because
%     % they are for different length.
%     [pos(1,1), scores(1,1)] = OVERLAPPING_DTW_MEX(theory,'experiment.txt', N, sets.comparison.R);
%     [pos(1,2), scores(1,2)] = OVERLAPPING_DTW_MEX(theory,'experimentrev.txt', N, sets.comparison.R);

%     % then we have one of the two choices for position
%     if scores(1)< scores(2)
%         rezMax.maxcoef = scores(1);
%         rezMax.pos = pos(1);
%         rezMax.or = 1;
%     else
%         rezMax.maxcoef = scores(2);
%         rezMax.pos = pos(2);
%         rezMax.or = 2;   
%     end
    % faster function, only when barC has bitmask only on left and right
    %     xcorrs = unmasked_pcc_corr(barC, theorBar, barB);
    %     [rezMax.maxcoef,rezMax.pos,rezMax.or] = get_best_parameters(xcorrs, 3 );
    %     % now find the maximum score for this stretching parameter
    %     xcorrMax(j) = rezMax{j}.maxcoef(1);
                    
end

%%
% 
%     % for query, only select the non-nan bit
%     query = query(logical(querybit));
%     locStart = find(query,1,'first');
%     
%     R= 5;  % describes the Sakoe-Chiba band, so in this case we don't allow any local stretching
%     % and instead increase this if we want better output
% %   
%     lS=1;
%     lD = 110000;
%     dTemp = data(lS:lD);
%     databitTem = databit(lS:lD);
% % mex 'NAN_DTW_MEX.cpp';
%     % check: what happens if bitmask is all ones?
% %     tic
%     [posD(1), score1] = NAN_DTW_MEX(dTemp,query,double(databitTem<50), length(dTemp), length(query), R);
%     [posD(2), score2] = NAN_DTW_MEX(dTemp,fliplr(query),double(dTemp<50), length(dTemp), length(query), R);
% %     toc
%     
%     % take the negative of minscore - then this will be max
%     [mincoef,or] = min([score1,score2]);
%     maxcoef = -mincoef;
%     pos = posD(or)+1; % if we want location for the full query
%     
%     secondPos = locStart;
%     lenM = length(query);
% 
%     %     
%     % relation of R to MAXSAMP?
%     %     seems ok!!
%     b1 = zscore(dTemp(pos:pos+length(query)-1),1);
%     b2 = zscore(query,1);
%     if or==2
%         b2=fliplr(b2);
%     end
% 
%     [val,IX,IY] = dtw(b1,b2,'squared',R);
%     
%     dist=[IX IY];
%     
%     mincoef
%     sqrt(val)
%     

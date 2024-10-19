function [sources] = track_sources(vec,distPar,maxDistFrames,Wsn)

% from hpfl
if nargin < 3
% not usaed really
maxDistFrames = inf;
distPar = 5;
end

if nargin < 4
   Wsn = vec;
end

% % source tracking (for regularization or multi-frame fluorophores
% vec = cell(1,100);
% for j=1:length(X)
%     for i=1:size(X{j},1)
%         vec{X{j}(i,1)} = [vec{X{j}(i,1)} X{j}(i,2)];
%     end
% end
% for i=1:size(X{2},1)
%     vec{X{2}(i,1)} = [vec{X{2}(i,1)} X{2}(i,2)];
% end

sources = cell(0,1);
init = 1;
while sum(vec{init})==0
    init=init+1;
end

for i = 1:length(vec{init})
    sources{length(sources)+1}=[init,vec{init}(i),Wsn{init}(i)]; % could use intensity
end



% distPar=5;
for l = init+1:length(vec)
    lastSources =cellfun(@(x) x(end,2),sources);
    lastSourcesIdx =cellfun(@(x) x(end,1),sources);
    if~isempty(vec{l})
        % maybe as first step use match-pairs between the previous time-frame,
        % and then pdist between the rest (possibly improves accuracy?)
        dists = (vec{l}'-lastSources').^2;
    %     A=pdist2(vec{l},lastSources'); % could be within some distance only allowed

        dists(l-lastSourcesIdx>maxDistFrames,:) = inf; % those that are larger than maxDist, don't assign
        dists(dists>distPar^2)=inf; %todo: adapt so it would be harder to ma over long gaps

            temp=zeros(length(vec{l}),1);

        M = matchpairs(dists,distPar^2/2);

        for a=1:size(M,1)
            sources{M(a,1)}=[sources{M(a,1)};[l,vec{l}(M(a,2)),Wsn{l}(M(a,2))]];
        end

        temp(M(:,2))=1;

        for i = 1:length(temp)
            if temp(i)==0
                sources{length(sources)+1}=[l,vec{l}(i),Wsn{l}(i)];
            end
        end
    end
end


% 
% for i=2:100
%     Xest =vec{i};
%     Xtrue =vec{2};
% 
%     dists = (Xest-Xtrue').^2;
%     M = matchpairs(dists, 2);
% 
% end
% 
%     M = matchpairs(dists, localWindowRadius);
%     
%     
% nblamb = length(lambda);
% sources = cell(0,1);
% 
% if nargin <4
%     distPar = 5; % distance so that features would be assigned to the same
% end
% 
% if nargin < 5
%     maxDistFrames = inf; % max distance between frames
% end
% 
% init = 1;
% while sum(X{init})==0
%     init=init+1;
% end
% 
% for i = 1:length(X{init})
%     sources{length(sources)+1}=[init,X{init}(i,:),amps{init}(i,:)];
% end
% 
% % distPar=5;
% for l = init+1:nblamb
%     lastSources =cellfun(@(x) x(end,2),sources);
%     lastSourcesIdx =cellfun(@(x) x(end,1),sources);
%     
%     % maybe as first step use match-pairs between the previous time-frame,
%     % and then pdist between the rest (possibly improves accuracy?)
% 
%     A=pdist2(X{l},lastSources'); % could be within some distance only allowed
%     
%     A(:,l-lastSourcesIdx>maxDistFrames) = inf; % those that are larger than maxDist, don't assign
%     A(A>distPar)=inf; %todo: adapt so it would be harder to ma over long gaps
% 
%     temp=zeros(length(X{l}(:,1)),1);
%     for a = 1:length(X{l}(:,1)) % go through all previous peaks
%         [val, n] = min(A(:));
%         if ~isinf(val)
%             [x,y] = ind2sub(size(A),n); % x-row/y-column
%             temp(x)=1; %this row has a match in l
%             % sources{x} might not hold the correct idx.. if there was some
%             % shift..
%             % check all last sources
% %             lastSources =cellfun(@(x) x(end,2),sources);
%             % compare to sources from previous step
% %             [val2,b] = min(pdist2(X{l}(x,:),lastSources')); % this gives an index in sources
%             %if val2 < distPar
% %                 xLast=b;
%                 % now update sources with the new element
%                 sources{y}=[sources{y};[lambda(l),X{l}(x,:),amps{l}(x,:)]];
%                 A(x,:)=Inf;
%                 A(:,y)=Inf;
% %             else
% %                 A(x,y)=Inf;
%         end
% 
%     end
% 
%     for i = 1:length(temp)
%         if temp(i)==0
%             sources{length(sources)+1}=[lambda(l),X{l}(i,:),amps{l}(i,:)];
%         end
%     end
% end

%
% 
% f=figure('Position', [100, 100, 600, 300])
% %  tiledlayout(2,1)
% %  nexttile
% % imagesc(kymo)
% %  nexttile
% hold on
% for s = 1:length(sources)
%     plot( sources{s}(:, 2),sources{s}(:, 1),'|','linewidth', 2)
% end
% % ylim([3.95, 5.05])
% set(gca, 'YDir','reverse')

end


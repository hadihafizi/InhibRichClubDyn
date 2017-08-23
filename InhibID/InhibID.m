
% function [FR, SpontFrac, InOutDeg, ISIstd] = InhibID(dataSet)
clear all
dataset = {'02-0','03-0','03-1','04-0','04-1','05-0','05-1',...
    '06-0','06-1','07-0','07-1','08-0','08-1','09-0','09-1'};
i = 1;
load([dataset{i},'/asdf.mat'])
TotT = asdf_raw{end}(2);
load([dataset{i},'/sig_del.mat'])
load([dataset{i},'/sig_te.mat'])
% [cShapes,cSizes,spontEvents] = asdf2avlets(asdf_raw,sig_del);
% save('spontEvents.mat','spontEvents')
load([dataset{i},'/spontEvents.mat'])
load('AllanFactors_Allinvitro.mat')

%% std(ISI) vs FR of spontaneous activity

SpontEventsISIStd = cell2mat(cellfun(@(x) std(diff(single(x))),spontEvents,'uni',false))/1000;
SpontEventsISI = (cellfun(@(x) diff(single(x))/1000,spontEvents,'uni',false));
SpontEventsISIAll = SpontEventsISI{1};
for i=2:length(SpontEventsISI)
    SpontEventsISIAll = [SpontEventsISIAll,SpontEventsISI{i}];
end
[ISInum,ISIbin] = hist(SpontEventsISIAll,50);
figure;loglog(ISIbin,ISInum)
SpontFR = cell2mat(cellfun(@length,spontEvents,'uni',false))/3600;

ISIstd = cell2mat(cellfun(@(x) std(diff(x)),asdf_raw(1:end-2),'uni',false))/1000;% std in seconds
ISI = (cellfun(@(x) diff(x)/1000,asdf_raw(1:end-2),'uni',false));% ISI in seconds
ISIAll = ISI{1};
for i=2:length(ISI)
    ISIAll = [ISIAll,ISI{i}];
end
[ISInum,ISIbin] = hist(ISIAll,50);
figure;loglog(ISIbin,ISInum)
FR = cell2mat(cellfun(@length,asdf_raw(1:end-2),'uni',false))/3600;

% figure; loglog(SpontFR,SpontEventsISIStd,'*','color','r'); hold on
% loglog(FR,ISIstd,'.','color','k')
% legend('Spontaneous','All')
% xlabel('Firing Rate [spike/s]')
% ylabel('std(ISI) [s]')


%% Neurons w/ High Spontaneous FR

[A1,B1] = sort(SpontFR,'descend'); num_hiSpont = ceil(0.2*length(A1));
HiSpontIdx = B1(1:num_hiSpont);  
%% Rich Neurons
% 
% load PDF_1_16_30ms.mat;load wgts_1_16ms.mat;
% sig_te = wgt.*PDF(:,:,45); oute = sum(sig_te,2); oute = oute./max(oute);
% [A1,B1] = sort(oute,'descend'); num_rich = ceil(0.2*length(A1));
% rich = B1(1:num_rich);  
% 
%% Correlation b/w Inhib and Rich
% 
% intersec = intersect(HiSpontIdx,rich);

%% std(ISI) vs FR vs (#indeg-#outdeg)/sum
% Hypothesis: inhibitory neurons have high firing rates, low variability of
% ISI, more incoming connection than outgoing ones (b/c TE cannot detect
% inhibitory connections) and high fraction of spontaneous activity (again
% measured by c-webs using connections defined by TE).

InOutDeg = (sum(sig_te,1) - sum(sig_te,2)')./ ...
    (sum(sig_te,1)+sum(sig_te,2)'); % difference b/w in- and out-deg of all neurons normalized by their sum
SpontFrac = SpontFR./FR; % fraction of spontaneous activity

figure;plot3k([ISIstd,SpontFrac,InOutDeg'],'Marker',{'.' 30});%grid on
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
% set(gca, 'ZScale', 'log')
% xlabel('indeg-outdeg')
% xlabel('Avg FR [spikes/s]')
ylabel('Fraction of Spontaneous Firing Rate')
zlabel('Indeg - Outdeg')
xlabel('std(ISI) [s]')

%% Clustering - All Data

ISI = cell(length(dataset),1); ISIstd = cell(length(dataset),1); 
ISI_All = [];
InOutDeg = cell(length(dataset),1);
FR = cell(length(dataset),1); SpontFR = cell(length(dataset),1);
SpontFrac = cell(length(dataset),1); AllPars = cell(length(dataset),1);
AllPars_Agg = [];

for i = 1:length(dataset)
    load([dataset{i},'/asdf.mat'])
    load([dataset{i},'/spontEvents.mat'])
    load([dataset{i},'/sig_te.mat'])
    ISI{i} = cellfun(@(x) diff(x)/1000,asdf_raw(1:end-2),'uni',false);% ISI in seconds
    ISIstd{i} = cell2mat(cellfun(@(x) std(diff(x)),asdf_raw(1:end-2),'uni',false))/1000;% std in seconds
    InOutDeg{i} = (sum(sig_te,1) - sum(sig_te,2)')./ ...
        (sum(sig_te,1)+sum(sig_te,2)'); % difference b/w in- and out-deg of all neurons normalized by their sum
    FR{i} = cell2mat(cellfun(@length,asdf_raw(1:end-2),'uni',false))/(asdf_raw{end}(2)/1000); % spks/s
    SpontFR{i} = cell2mat(cellfun(@length,spontEvents,'uni',false))/(asdf_raw{end}(2)/1000);
    SpontFrac{i} = SpontFR{i}./FR{i}; % fraction of spontaneous activity
    
    ISI_All = [ISI_All; ISI{i}];
    AF{i}(isnan(AF{i})) = 0;
%     AllPars{i} = [AF{i}',ISIstd{i},SpontFR{i},SpontFrac{i},InOutDeg{i}'];
%     AllPars{i} = [AF{i}',ISIstd{i},SpontFR{i},SpontFrac{i}];
    AllPars{i} = [AF{i}(26,:)',ISIstd{i}];
%     AllPars{i} = [AF{i}'];
    AllPars_Agg = [AllPars_Agg;AllPars{i}]; % Aggregate All parameters from all data sets
end

opts = statset('Display','final');
[idx,C] = kmeans(AllPars_Agg,2,'Distance','cityblock',...
    'Replicates',5,'Options',opts);

figure;
plot3(AllPars_Agg(idx==1,end-3),AllPars_Agg(idx==1,end-2),AllPars_Agg(idx==1,end),'r.','MarkerSize',20)
hold on
plot3(AllPars_Agg(idx==2,end-3),AllPars_Agg(idx==2,end-2),AllPars_Agg(idx==2,end),'b.','MarkerSize',20)
plot3(AllPars_Agg(idx==3,end-3),AllPars_Agg(idx==3,end-2),AllPars_Agg(idx==3,end),'k.','MarkerSize',20)
xlabel 'std(ISI) [s]'
ylabel 'Spontaneous FR [spk/s]'
zlabel '(InDeg-OutDeg)/sum'
legend('Cluster 1','Cluster 2','Location','NW')
title 'Cluster Assignments'
hold off
grid on

%%% plot ISI Distributions for each category
idx1 = 1:size(AllPars_Agg,1); idx2 = 1:size(AllPars_Agg,1);
idx1 = idx1(idx==1); idx2 = idx2(idx==2);
ISI_Agg1 = []; ISI_Agg2 = [];
for k = 1:length(idx1)
    ISI_Agg1 = [ISI_Agg1, ISI_All{idx1(k)}];
end
for k = 1:length(idx2)
    ISI_Agg2 = [ISI_Agg2, ISI_All{idx2(k)}];
end
figure;
[idx1counts, idx1centers] = hist(ISI_Agg1,(0:max(ISI_Agg1))); idx1counts = idx1counts/sum(idx1counts);
[idx2counts, idx2centers] = hist(ISI_Agg2,(0:max(ISI_Agg2))); idx2counts = idx2counts/sum(idx2counts);
semilogx(idx1centers, idx1counts, 'r-*', 'LineWidth', 2); hold on
semilogx(idx2centers, idx2counts, 'b-*', 'LineWidth', 2);
legend('idx1','idx2')
title 'ISI Distributions for K-means Categories'
ylabel 'Probability'
xlabel 'ISI [s]'

%%%%% PCA
% AllParstmp = [AF{1}(26,:)',ISIstd,SpontFrac,InOutDeg'];
% [coeff,score,latent,tsquared,explained,mu] = pca(AllParstmp);
% 
% figure;
% plot(score(:,1),score(:,2),'b.','MarkerSize',20);


%% t-SNE on All Parameters


tsneSpace = tsne(AllPars_Agg);

figure;
plot(tsneSpace(idx==1,1),tsneSpace(idx==1,2),'b.','MarkerSize',10); hold on
plot(tsneSpace(idx==2,1),tsneSpace(idx==2,2),'r.','MarkerSize',10)
plot(tsneSpace(idx==3,1),tsneSpace(idx==3,2),'k.','MarkerSize',10)
title 't-SNE Space- AF(6s), std(ISI)'

%% t-SNE on ISI Distributions

ISIcounts = zeros(length(ISI_All), length(0:max(cellfun(@max, ISI_All))));
ISIcenters = zeros(1, length(0:max(cellfun(@max, ISI_All))));
for k = 1:length(ISI_All)
    [ISIcounts(k,:), ISIcenters] = hist(ISI_All{k},(0:max(cellfun(@max, ISI_All)))); 
    ISIcounts(k,:) = ISIcounts(k,:)/sum(ISIcounts(k,:));
end

ISItsneSpace = tsne(ISIcounts);

[idx,C] = kmeans(ISIcounts,2,'Distance','cityblock',...
    'Replicates',5,'Options',opts);

figure;
plot(ISItsneSpace(idx==1,1),ISItsneSpace(idx==1,2),'b.','MarkerSize',10); hold on
plot(ISItsneSpace(idx==2,1),ISItsneSpace(idx==2,2),'r.','MarkerSize',10)
plot(ISItsneSpace(idx==3,1),ISItsneSpace(idx==3,2),'k.','MarkerSize',10)
title 't-SNE Space- ISI Dist'

figure;
semilogx(ISIcenters(1,:), mean(ISIcounts(idx==1,:)), 'r-*', 'LineWidth', 2); hold on
semilogx(ISIcenters(1,:), mean(ISIcounts(idx==2,:)), 'b-*', 'LineWidth', 2);
legend('idx1','idx2')
title 'ISI Distributions for K-means Categories'
ylabel 'Probability'
xlabel 'ISI [s]'

%% t-SNE on JSDivergences of ISI Distributions

asdf_all = cell(0,0);
for i = 1:length(dataset)
    load([dataset{i},'/asdf.mat']);
    asdf_all = [asdf_all;asdf_raw(1:end-2)];
end
asdf_all = [asdf_all;asdf_raw(end-1:end)];
asdf_all{end}(1) = length(asdf_all) - 2;

[ jsdMat ] = ssJSD_DistMatrix_Hadi ( asdf_all, [.1 .5 1 .5 .1] );
jsdMat = jsdMat + jsdMat';

tsneSpace = tsne(jsdMat, 3);
idx = kmeans(jsdMat,2);

figure;
plot(tsneSpace(idx==1,1),tsneSpace(idx==1,2),'b.','MarkerSize',10); hold on
plot(tsneSpace(idx==2,1),tsneSpace(idx==2,2),'r.','MarkerSize',10);
plot(tsneSpace(idx==3,1),tsneSpace(idx==3,2),'k.','MarkerSize',10);
plot(tsneSpace(idx==4,1),tsneSpace(idx==4,2),'y.','MarkerSize',10);
plot(tsneSpace(idx==5,1),tsneSpace(idx==5,2),'c.','MarkerSize',10);



%% Correlating different parameters
% 
% topPercent = 0.35; % the percentage of the data used.
% [A1,B1] = sort(oute,'descend'); num_rich = ceil(topPercent*length(A1));
% rich = B1(1:num_rich);  
% [A1,B1] = sort(FR,'descend'); num_hiFR = ceil(topPercent*length(A1));
% HiFRIdx = B1(1:num_hiFR);  
% [A1,B1] = sort(SpontFR,'descend'); num_hiSpont = ceil(topPercent*length(A1));
% HiSpontIdx = B1(1:num_hiSpont);  
% [A1,B1] = sort(ISIstd,'ascend'); num_loISI = ceil(topPercent*length(A1));
% loISIIdx = B1(1:num_loISI);  
% [A1,B1] = sort(InOutDeg,'descend'); num_hiInOutDeg = ceil(topPercent*length(A1));
% HiInOutIdx = B1(1:num_hiInOutDeg);  
% 
% FRnISI = intersect(HiFRIdx,loISIIdx);
% FRnInOut = intersect(HiFRIdx,HiInOutIdx);
% FRnSpontFrac = intersect(HiFRIdx,HiSpontIdx);
% 
% ISInSpontFrac = intersect(loISIIdx,HiSpontIdx);
% ISInInOut = intersect(loISIIdx,HiInOutIdx);
% 
% InOutnSpontFrac = intersect(HiInOutIdx,HiSpontIdx);
% 
% MajorIntersects = intersect(FRnISI,...
%     intersect(FRnSpontFrac,ISInSpontFrac));
%     
% AllIntersect = intersect(FRnISI,...
%     intersect(FRnInOut,...
%     intersect(FRnSpontFrac,...
%     intersect(ISInSpontFrac,...
%     intersect(ISInInOut,InOutnSpontFrac)))));



clear all

suffs = {'rand','vweakRC','weakRC','medRC1','medRC2','strongRC','rwNoRC'};
% parpool(10);
%% RC coefficient

for i = 1:length(suffs)
    load(['Simk100_Generic_SamePoissSeeds_rcWtNetTotInOutMethod2_',suffs{i},'_500nodes.mat'], ...
        'weightmat')
    Rw = rich_club_wd(weightmat);
    RwR = zeros(500,length(Rw));
    parfor k = 1:500
        R = randmio_dir(weightmat,10);
        RwR(k,:) = rich_club_wd(R);
    end
    
    NodeDegree = sum((weightmat~=0))+sum((weightmat'~=0));
    figure(104);
    plot(1:max(NodeDegree), Rw./mean(RwR,1),'LineWidth',2);hold on
    
end
xlabel('Node Degree (K)')
ylabel('Normalized Weighted RC Coef')
legend(suffs)

%% Firing rate and Avalanche length/size distributions

numAvals = zeros(length(suffs), 1);
AvgFR = zeros(length(suffs), 1); stdFR = zeros(length(suffs), 1);
AvgAvalLen = zeros(length(suffs), 1); AvgAvalSize = zeros(length(suffs), 1);
stdAvalLen = zeros(length(suffs), 1); stdAvalSize = zeros(length(suffs), 1);

for i = 1:length(suffs)
    load(['Simk100_Generic_SamePoissSeeds_rcWtNetTotInOutMethod2_',suffs{i},'_500nodes.mat'], ...
        'asdf')
    
    firingRate = cellfun(@numel,asdf_raw(1:asdf_raw{end}(1)))/(asdf_raw{end}(2)/1000); % spikes/s
    [FRcounts, FRedges] = histcounts(firingRate);
    FRcounts = FRcounts/(sum(FRcounts)*(FRedges(2) - FRedges(1)));
    FRcenters = movsum(FRedges, 2)/2;
    FRcenters = FRcenters(2:end);
    figure(105);
    loglog(FRcenters, FRcounts, 'LineWidth', 2); hold on
    xlabel('firing rate [spike/s]')
    ylabel('Probability Density')
    
    [aShapes,aLengths,aSizes] = asdf2Shapes(asdf);
    [Lencounts,Lenedges] = histcounts(aLengths);
    Lencounts = Lencounts/(sum(Lencounts)*(Lenedges(2) - Lenedges(1)));
    Lencenters = movsum(Lenedges, 2)/2;
    Lencenters = Lencenters(2:end);
    figure(106);
    loglog(Lencenters, Lencounts,'LineWidth',2);hold on
    xlabel('Avalanche Length [ms]')
    ylabel('Probability Density')
    
    [Sizecounts,Sizeedges] = histcounts(aSizes);
    Sizecounts = Sizecounts/(sum(Sizecounts)*(Sizeedges(2) - Sizeedges(1)));
    Sizecenters = movsum(Sizeedges, 2)/2;
    Sizecenters = Sizecenters(2:end);
    figure(107);
    loglog(Sizecenters, Sizecounts,'LineWidth',2);hold on
    xlabel('Avalanche Sizes [#neurons]')
    ylabel('Probability Density')
    
    numAvals(i) = numel(aLengths);
    AvgFR(i) = mean(firingRate); stdFR(i) = std(firingRate);
    AvgAvalLen(i) = mean(aLengths); stdAvalLen(i) = std(aLengths);
    AvgAvalSize(i) = mean(aSizes); stdAvalSize(i) = std(aSizes);
    
end

figure(105); legend(suffs)
figure(106); legend(suffs)
figure(107); legend(suffs)

figure;
bar(AvgFR, 'LineWidth', 2); hold on
errorbar(AvgFR, stdFR, 'kx', 'LineWidth', 2)
title('Firing Rate')
ax = gca; ax.XTick = 1:length(AvgFR);
ax.XTickLabel = suffs;

figure;
bar(numAvals, 'LineWidth', 2); hold on
title('Number of Avalanches')
ax = gca; ax.XTick = 1:length(AvgFR);
ax.XTickLabel = suffs;

figure;
bar(AvgAvalLen, 'LineWidth', 2); hold on
errorbar(AvgAvalLen, stdAvalLen, 'kx', 'LineWidth', 2)
title('Avg Avalanche Length')
ax = gca; ax.XTick = 1:length(AvgFR);
ax.XTickLabel = suffs;

figure;
bar(AvgAvalSize, 'LineWidth', 2); hold on
errorbar(AvgAvalSize, stdAvalSize, 'kx', 'LineWidth', 2)
title('Avg Avalanche Size')
ax = gca; ax.XTick = 1:length(AvgFR);
ax.XTickLabel = suffs;








% avalanche synchrony
% calculate Jaccard coefficients for neurons participating in avalanches of
% different lengths

clear all

dataSet = {'02-0','03-0','03-1','04-0','04-1','05-0','05-1','06-0',...
    '06-1','07-0','07-0','08-0','08-1','09-0','09-1'};
% dataSet = {'139', '151', '152', '163', '165', '168', '174'};
% dataSet = {'rand','vweakRC','weakRC','medRC1','medRC2','strongRC','rwNoRC'};
PFeig = [102,102,103,104,104,107,103];
BinSize = 1;
len = [2 5 8 12 15 20 25 30];

%%
for i = 1:length(dataSet)
    tic
    asdf = load([dataSet{i},'/asdf.mat']);
%     asdf = load([dataSet{i},'/asdf_rest',dataSet{i},'_TimeCorOneMSGap.mat']);
    asdf.asdf_raw = ASDFChangeBinning(asdf.asdf_raw, BinSize);
%     asdf = load(['Simk',num2str(PFeig(i)),'_Generic_',dataSet{i},'vLowerSeeds.mat'],'asdf');
    str_aval = asdf2Shapes(asdf.asdf_raw);
    sim = zeros(length(asdf.asdf_raw)-2, length(asdf.asdf_raw)-2, length(len));
    for ilen = 1:length(len)
        avals_id = cellfun(@(x) length(unique(x(:,2))) == len(ilen), str_aval);
        avals_len = str_aval(avals_id);
        new_asdf = cell(length(asdf.asdf_raw)-2, 1);
        for anum = 1:length(avals_len)
            for neu=1:length(avals_len{anum})
                new_asdf{avals_len{anum}(neu,1)} = [new_asdf{avals_len{anum}(neu,1)}, avals_len{anum}(neu,2)];
            end
        end
        new_asdf{end+1} = 1;
        new_asdf{end+1} = [length(asdf.asdf_raw)-2, 3.6e6];
        sim(:,:,ilen) = sim_spktrn(new_asdf);
    end
%     parsave([dataSet{i},'\AvalSynchrony_',dataSet{i},'_bs',num2str(BinSize),'ms.mat'],{'sim'},{sim})
%     parsave(['AvalSynchrony_',dataSet{i},'.mat'],{'sim'},{sim})

    toc
end

%% plot synchrony distributions

sim_len = cell(length(len),1);
for i = 1:length(dataSet)
    load([dataSet{i},'/AvalSynchrony_',dataSet{i},'_bs',num2str(BinSize),'ms.mat'])
%     load(['AvalSynchrony_',dataSet{i},'.mat'])
    for ilen = 1:length(len)
        temp = sim(:,:,ilen);
        sim_len{ilen} = [sim_len{ilen}; temp(triu(temp, 1)~=0)];
    end
end


figure;
subplot(2,2,1);
for ilen = 1:length(len)

    [~, edges] = histcounts(sim_len{ilen}, 'Normalization', 'probability');
    [rrcount, ~] = histcounts(sim_len{ilen}, edges, 'Normalization', 'probability');
    centers = movsum(edges, 2)/2;
    centers = centers(2:end);

% plot(centers,rrcount,'r--.','MarkerSize',10); hold on
    semilogx(centers,rrcount,'LineWidth',1); hold on
end
xlabel('Jaccard Coef')
ylabel('Probablity')
legend('2','5','8','12','15','20','25','30')
set(gca,'tickDir','out','color','none','fontname','arial','linewidth',1,'fontsize',9); box off;

for ilen = 1:length(len)
    sim_len{ilen}(isnan(sim_len{ilen})) = [];
end
avg_sim = cellfun(@mean, sim_len);
median_sim = cellfun(@median, sim_len);
std_sim = cellfun(@std, sim_len);
subplot(2,2,2);
f1 = shadedErrorBar(len, avg_sim, std_sim/sqrt(length(dataSet))); hold on
% f2 = plot(len, median_sim, 'LineWidth', 1);
xlabel('Avalanche Length [ms]')
ylabel('Jaccard Coef')
% legend([f1.mainLine, f2],'Mean', 'Median', 'Location', 'Best')
set(gca,'tickDir','out','color','none','fontname','arial','linewidth',1,'fontsize',9); box off;

%% corr

min_numel = numel(sim_len{end});
% r = zeros(length(len)-1, 1);
% pvalue = zeros(length(len)-1, 1);
sim_len_mat = [];
for ilen = 1:length(len)
    sim_len_mat = [sim_len_mat, sim_len{ilen}(randsample(length(sim_len{ilen}), min_numel))];
end
[r, pvalue] = corrcoef(sim_len_mat);
%% ANOVA
varname = cell(min_numel, 1);
varname(:) = {'sim'};
t = table(varname, sim_len_mat(:,1), sim_len_mat(:,2), sim_len_mat(:,3), sim_len_mat(:,4), ...
    sim_len_mat(:,5), sim_len_mat(:,6), sim_len_mat(:,7), sim_len_mat(:,8), ...
    'VariableNames',{'AvalLen','len2','len5','len8','len12','len15','len20','len25','len30'});
Meas = dataset([1:8]','VarNames',{'Measurements'});
rm = fitrm(t,'len2-len30~AvalLen','WithinDesign',Meas);
anova(rm)

%% t-test

[h, p] = ttest(sim_len_mat(:,1), sim_len_mat(:,2), 'Tail', 'right')



%%

% r2 = cellfun(@(x,y) xcorr(x,y,0), sim_len, sim_len);

figure;
plot(len(1:end-1), r); hold on
% plot(len, r2);
set(gca,'tickDir','out','color','none','fontname','arial','linewidth',1,'fontsize',9); box off;




































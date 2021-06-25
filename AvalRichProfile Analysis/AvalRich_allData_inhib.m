

% function AvalRichProfileRev_DoItAll(dataset,BinSize)

% This script loads actual and shuffled asdf files for the given data
% set and generates avalanche string files "aval_listXXX_bsXms.mat" for all
% of them and also avalanche richness profiles "AvalRichProfile_bsXms.mat"
% and then uses those data to plot the avalanche richness profile of actual
% and averaged shuffled data sets for multiple avalanche lengths.
% 
% Requires functions: Avalanche_detection.m
%                     AvalRichProfileRev.m
%                     AvgCollapse.m
%                     shadedErrorBar.m
% 
% Requires files: asdf.mat
%                 PDF_1_16_30ms.mat or PDF_NoSpur.mat
%                 wgts_1_16ms.mat
%                 aval_listXXX.mat (if using saved files)
% 
% Hadi Hafizi, Jan. 2016



clear all
% dataset = {'02-0','03-0','03-1','04-0','04-1','05-0','05-1','06-0',...
%     '06-1','07-0','08-0','08-1','09-0','09-1'};
dataset = {'139', '151', '152', '163', '165', '168', '174'};
EIdatasets = [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25];
suffix = '_noUnID_cor_toteRC';
BinSize = 3; %ms

len = [2 5 8 12 15 20 25 30];
% 
rchnss_all = {[],[],[],[],[],[],[],[],[],[],[],[],[],[]};
err_all = {[],[],[],[],[],[],[],[],[],[],[],[],[],[]};
rchnss_nr_all = {[],[],[],[],[],[],[],[],[],[],[],[],[],[]};
err_nr_all = {[],[],[],[],[],[],[],[],[],[],[],[],[],[]};
rchnss_nre_all = {[],[],[],[],[],[],[],[],[],[],[],[],[],[]};
err_nre_all = {[],[],[],[],[],[],[],[],[],[],[],[],[],[]};
rchnss_nri_all = {[],[],[],[],[],[],[],[],[],[],[],[],[],[]};
err_nri_all = {[],[],[],[],[],[],[],[],[],[],[],[],[],[]};
rchnss_re_all = {[],[],[],[],[],[],[],[],[],[],[],[],[],[]};
err_re_all = {[],[],[],[],[],[],[],[],[],[],[],[],[],[]};
rchnss_ri_all = {[],[],[],[],[],[],[],[],[],[],[],[],[],[]};
err_ri_all = {[],[],[],[],[],[],[],[],[],[],[],[],[],[]};
rchnss_mean = {[],[],[],[],[],[],[],[]};
err_mean = {[],[],[],[],[],[],[],[]};
rchnss_nr_mean = {[],[],[],[],[],[],[],[]};
err_nr_mean = {[],[],[],[],[],[],[],[]};
rchnss_nre_mean = {[],[],[],[],[],[],[],[]};
err_nre_mean = {[],[],[],[],[],[],[],[]};
rchnss_nri_mean = {[],[],[],[],[],[],[],[]};
err_nri_mean = {[],[],[],[],[],[],[],[]};
rchnss_re_mean = {[],[],[],[],[],[],[],[]};
err_re_mean = {[],[],[],[],[],[],[],[]};
rchnss_ri_mean = {[],[],[],[],[],[],[],[]};
err_ri_mean = {[],[],[],[],[],[],[],[]};

numR = zeros(length(dataset), 1);
numNR = zeros(length(dataset), 1);
numR_E = zeros(length(dataset), 1);
numR_I = zeros(length(dataset), 1);
numNR_E = zeros(length(dataset), 1);
numNR_I = zeros(length(dataset), 1);

for i = 1:length(dataset)
    %
    PDF = load([dataset{i},'/PDF_NoSpur_thr45_',dataset{i},'.mat']);
    wgt = load([dataset{i},'/wgts_1_16ms.mat']);
    sig_te = wgt.wgt.*PDF.PDF_cor;
    
    exc_thresh = .9;
    inh_thresh = .9;
    % in-vitro
%     neurontype = load(['C:\Users\Mahzad\Box Sync\Research\Beggs Lab\Projects\InhibID\Ian Stevenson lab\Neuron Type Analysis\Neuron Type Analysis\in vitro\DataSet', num2str(EIdatasets(i)), '.mat']);
    % in-vivo
    neurontype = load(['C:\Users\Mahzad\Box Sync\Research\Beggs Lab\Projects\InhibID\Ian Stevenson lab\Neuron Type Analysis\Neuron Type Analysis\in vivo\unitactivity_Mouse ', num2str(dataset{i}), '.mat']);
    excprobs = extractfield(neurontype.neurontype, 'excprob');
    inhprobs = extractfield(neurontype.neurontype,'inhprob');
    % excConRatio = extractfield(neurontype, 'ExcConRatio');
    ExcIdx = find(excprobs > exc_thresh);
    InhIdx = find(inhprobs > inh_thresh);
    neuron_ids = 1:numel(neurontype.neurontype);
    unidentified = neuron_ids(~ismember(neuron_ids, ExcIdx) & ~ismember(neuron_ids, InhIdx));
    
    % filter out the Unidentified neurons
    excprobs = excprobs(~ismember(1:length(excprobs), unidentified));
    inhprobs = inhprobs(~ismember(1:length(inhprobs), unidentified));
    ExcIdx = find(excprobs > exc_thresh);
    InhIdx = find(inhprobs > inh_thresh);

    sig_te = sig_te(~ismember(1:length(sig_te), unidentified), ~ismember(1:length(sig_te), unidentified));
    oute = sum(sig_te,2);
    inte = sum(sig_te,1);
    tote = oute + inte';
    [A1,B1] = sort(tote,'descend');
    num_rich = ceil(0.2*length(A1));
    r_idx = B1(1:num_rich);
    nr_idx = B1(num_rich+1:end);
%
    numR(i) = 1;%numel(r_idx);
    numNR(i) = 1;%numel(nr_idx);
    numR_E(i) = 1;%numel(intersect(r_idx, ExcIdx, 'stable'));
    numR_I(i) = 1;%numel(intersect(r_idx, InhIdx, 'stable'));
    numNR_E(i) = 1;%numel(intersect(nr_idx, ExcIdx, 'stable'));
    numNR_I(i) = 1;%numel(intersect(nr_idx, InhIdx, 'stable'));

    load([dataset{i},'/AvalRichProfile_',dataset{i},suffix,'_bs',num2str(BinSize),'ms.mat']);
    load([dataset{i},'/AvalRichProfile_',dataset{i},suffix,'_bs',num2str(BinSize),'ms_NR.mat']);
    load([dataset{i},'/AvalRichProfile_',dataset{i},suffix,'_bs',num2str(BinSize),'ms_NRE.mat']);
    load([dataset{i},'/AvalRichProfile_',dataset{i},suffix,'_bs',num2str(BinSize),'ms_NRI.mat']);
    load([dataset{i},'/AvalRichProfile_',dataset{i},suffix,'_bs',num2str(BinSize),'ms_RE.mat']);
    load([dataset{i},'/AvalRichProfile_',dataset{i},suffix,'_bs',num2str(BinSize),'ms_RI.mat']);
    AvalLen_all = AvalLen;
    rchnss_all{i} = cellfun(@(x) x./numR(i), rchnss, 'UniformOutput', false);
    err_all{i} = cellfun(@(x) x./numR(i), err, 'UniformOutput', false);
    rchnss_nr_all{i} = cellfun(@(x) x./numNR(i), rchnss_nr, 'UniformOutput', false);
    err_nr_all{i} = cellfun(@(x) x./numNR(i), err_nr, 'UniformOutput', false);
    rchnss_nre_all{i} = cellfun(@(x) x./numNR_E(i), rchnss_nre, 'UniformOutput', false);
    err_nre_all{i} = cellfun(@(x) x./numNR_E(i), err_nre, 'UniformOutput', false);
    rchnss_nri_all{i} = cellfun(@(x) x./numNR_I(i), rchnss_nri, 'UniformOutput', false);
    err_nri_all{i} = cellfun(@(x) x./numNR_I(i), err_nri, 'UniformOutput', false);
    rchnss_re_all{i} = cellfun(@(x) x./numR_E(i), rchnss_re, 'UniformOutput', false);
    err_re_all{i} = cellfun(@(x) x./numR_E(i), err_re, 'UniformOutput', false);
    rchnss_ri_all{i} = cellfun(@(x) x./numR_I(i), rchnss_ri, 'UniformOutput', false);
    err_ri_all{i} = cellfun(@(x) x./numR_I(i), err_ri, 'UniformOutput', false);
end
for ilen = 1:length(len)
    rchnss_mean{ilen} = zeros(size(rchnss_all{1}{ilen}));
    err_mean{ilen} = zeros(size(err_all{1}{ilen}));
    rchnss_nr_mean{ilen} = zeros(size(rchnss_nr_all{1}{ilen}));
    err_nr_mean{ilen} = zeros(size(err_nr_all{1}{ilen}));
    rchnss_nre_mean{ilen} = zeros(size(rchnss_nre_all{1}{ilen}));
    err_nre_mean{ilen} = zeros(size(err_nre_all{1}{ilen}));
    rchnss_nri_mean{ilen} = zeros(size(rchnss_nri_all{1}{ilen}));
    err_nri_mean{ilen} = zeros(size(err_nri_all{1}{ilen}));
    rchnss_re_mean{ilen} = zeros(size(rchnss_re_all{1}{ilen}));
    err_re_mean{ilen} = zeros(size(err_re_all{1}{ilen}));
    rchnss_ri_mean{ilen} = zeros(size(rchnss_ri_all{1}{ilen}));
    err_ri_mean{ilen} = zeros(size(err_ri_all{1}{ilen}));
end
for i = 1:length(dataset)
    for ilen = 1:length(len)
        rchnss_mean{ilen} = rchnss_mean{ilen} + rchnss_all{i}{ilen}/length(dataset);
        err_mean{ilen} = err_mean{ilen} + err_all{i}{ilen}/length(dataset);
        rchnss_nr_mean{ilen} = rchnss_nr_mean{ilen} + rchnss_nr_all{i}{ilen}/length(dataset);
        err_nr_mean{ilen} = err_nr_mean{ilen} + err_nr_all{i}{ilen}/length(dataset);
        rchnss_nre_mean{ilen} = rchnss_nre_mean{ilen} + rchnss_nre_all{i}{ilen}/length(dataset);
        err_nre_mean{ilen} = err_nre_mean{ilen} + err_nre_all{i}{ilen}/length(dataset);
        rchnss_nri_mean{ilen} = rchnss_nri_mean{ilen} + rchnss_nri_all{i}{ilen}/length(dataset);
        err_nri_mean{ilen} = err_nri_mean{ilen} + err_nri_all{i}{ilen}/length(dataset);
        rchnss_re_mean{ilen} = rchnss_re_mean{ilen} + rchnss_re_all{i}{ilen}/length(dataset);
        err_re_mean{ilen} = err_re_mean{ilen} + err_re_all{i}{ilen}/length(dataset);
        rchnss_ri_mean{ilen} = rchnss_ri_mean{ilen} + rchnss_ri_all{i}{ilen}/length(dataset);
        err_ri_mean{ilen} = err_ri_mean{ilen} + err_ri_all{i}{ilen}/length(dataset);
    end
end
%%
% 
% plot the figures for all the chosen lengths
rchnss = {rchnss_re_mean, rchnss_ri_mean, rchnss_nre_mean, rchnss_nri_mean};
err = {err_mean, err_re_mean, err_ri_mean, err_nr_mean, err_nre_mean, err_nri_mean};
avg_rchnss = zeros(length(rchnss), length(len));
color = {'r.', 'g', 'm', 'b', 'c', 'k','r', 'y'};
titles = {'Excitatory RC Participation','Inhibitory RC Participation',...
    'Excitatory NonRC Participation','Inhibitory NonRC Participation'};
H = figure;clf;
for i = 1:length(rchnss)
    subplot(2,2,i);
    for ilen = 2 : length(len) - 2
        %     set(H,'units','normalized','outerposition',[0 0 1 1])
        %     title('Proportion of Active Rich Neurons for Actual (Red) and Shuffled (Black) Data')
        
        f1 = shadedErrorBar(AvalLen{ilen},rchnss{i}{ilen},err{i}{ilen}/sqrt(length(dataset)), color{ilen});hold on
        avg_rchnss(i, ilen) = mean(rchnss{i}{ilen});
        %     f2 = shadedErrorBar(AvalLen{ilen},rchnss_re_mean{ilen},err_re_mean{ilen},'g');hold on
        %     f3 = shadedErrorBar(AvalLen{ilen},rchnss_ri_mean{ilen},err_ri_mean{ilen},'m');hold on
        %     f4 = shadedErrorBar(AvalLen{ilen},rchnss_nr_mean{ilen},err_nr_mean{ilen},'b');hold on
        %     f5 = shadedErrorBar(AvalLen{ilen},rchnss_nre_mean{ilen},err_nre_mean{ilen},'c');hold on
        %     f6 = shadedErrorBar(AvalLen{ilen},rchnss_nri_mean{ilen},err_nri_mean{ilen},'k');hold on
        %     f7 = plot(AvalLen{ilen},.3*ones(1,length(AvalLen{ilen})),'y');

    end
    grid on
    axis([0 len(6) 0 .4])
    title(titles{i})
    xlabel(['avalanche length [ms]']);
    % %     ylabel('mean # of active rich nodes')
    set(gca,'tickDir','out','color','none','fontname','arial','linewidth',1,'fontsize',12); box off;
%     xticks([3.33,6.66,10,13.33,16.66,20])
%     xticks([3.33,6.66,10,13.33,16.66,20])
%     xticklabels({'10','20','30','40','50','60'})
    xticklabels({'0','15','30','45','60'})
end
% legend([f1.mainLine, f2.mainLine, f3.mainLine, f4.mainLine, f5.mainLine, f6.mainLine],'R', 'RE', 'RI', 'NR', 'NRE', 'NRI')

% % save([dataset,'\AvalRichProfileShuf_AllDatasets_bs',num2str(BinSize),'ms.mat'],'AvalancheLen','RichnessShuff','ErrorShuff')
% % save([dataset,'\AvalRichProfile_AllDatasets_bs',num2str(BinSize),'ms.mat'],'AvalancheLen','Richness','Error')
% savefig([dataset,'\AvalRichActShuff_AllDatasets_bs',num2str(BinSize),'ms_subvivo'])
% saveas(H,[dataset,'\AvalRichActShuff_AllDatasets_bs',num2str(BinSize),'ms_subvivo.jpg'])

%%

R_IoverE = avg_rchnss(2,:)./avg_rchnss(1,:);
NR_IoverE = avg_rchnss(4,:)./avg_rchnss(3,:);
figure(30);clf;
% subplot(2,2,1);
plot(len(1:end), R_IoverE,'linewidth',1); hold on
plot(len(1:end), NR_IoverE,'linewidth',1);
% xticklabels({'15','30','45','60','75','90'})
xlabel('avalanche length [ms]')
ylabel({['I/E average avalanche'], ['participation ratio']})
legend('RC', 'nonRC')
set(gca,'tickDir','out','color','none','fontname','arial','linewidth',1,'fontsize',12); box off;

%% Statistical Test

re = cell(length(len), 1);
ri = cell(length(len), 1);
nre = cell(length(len), 1);
nri = cell(length(len), 1);

for ilen = 1:length(len)
    for i = 1:length(dataset)
        re{ilen} = [re{ilen}, rchnss_re_all{i}{ilen}];
        ri{ilen} = [ri{ilen}, rchnss_ri_all{i}{ilen}];
        nre{ilen} = [nre{ilen}, rchnss_nre_all{i}{ilen}];
        nri{ilen} = [nri{ilen}, rchnss_nri_all{i}{ilen}];
    end
end

% KS tests between consecutive avalanche lengths
re_test = zeros(length(len) - 1, 2);
ri_test = zeros(length(len) - 1, 2);
nre_test = zeros(length(len) - 1, 2);
nri_test = zeros(length(len) - 1, 2);
for ilen = 1:length(len) - 1
    [re_test(ilen, 1), re_test(ilen, 2)] = kstest2(re{ilen}, re{ilen + 1});
    [ri_test(ilen, 1), ri_test(ilen, 2)] = kstest2(ri{ilen}, ri{ilen + 1});
    [nre_test(ilen, 1), nre_test(ilen, 2)] = kstest2(nre{ilen}, nre{ilen + 1});
    [nri_test(ilen, 1), nri_test(ilen, 2)] = kstest2(nri{ilen}, nri{ilen + 1});
end




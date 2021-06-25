

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
dataset = '03-0';%[139 151 152 163 165 168 174];

BinSize = 1; %ms
load([dataset,'/asdf.mat'])

str_aval = Avalanche_detection(BinSize,asdf);
thr = 45; % PDF significance level (threshold)
% load([dataset,'\PDF_1_16_30ms.mat']);
% sig_te = wgt.*PDF(:,:,thr); 
load([dataset,'/PDF_NoSpur_thr45_',dataset,'.mat'])
load([dataset,'/wgts_1_16ms.mat']);
sig_te = wgt.*PDF; 
len = [2 5 8 12 15 20 25 30];
[AvalLen,rchnss,err] = AvalRichProfileRev(len,str_aval,sig_te,0);
% 
parfor shufnum = 1:20
%     x = load([dataset,'\asdf',num2str(shufnum),'.mat']);
%     str_aval = Avalanche_detection(BinSize,x.asdf);
    x = load([dataset,'/aval_list',dataset,'shuf',num2str(shufnum),'_bs',num2str(BinSize),'ms.mat']);
    [AvalLen_shuf,rchnss_shuf{shufnum},err_shuf{shufnum}] = ...
        AvalRichProfileRev(len,x.str_aval,sig_te);
end
for ilen = 1:length(len)
    rchnss_shuf_mean{ilen} = zeros(size(rchnss_shuf{1}{ilen}));
    err_shuf_mean{ilen} = zeros(size(err_shuf{1}{ilen}));
end
for shufnum = 1:20
    for ilen = 1:length(len)
        rchnss_shuf_mean{ilen} = rchnss_shuf_mean{ilen} + rchnss_shuf{shufnum}{ilen}/20;
        err_shuf_mean{ilen} = err_shuf_mean{ilen} + err_shuf{shufnum}{ilen}/20;
    end
end
%%
% 
% plot the figures for all the chosen lengths
H = figure(304);clf
for ilen = 1 : length(len)
    set(H,'units','normalized','outerposition',[0 0 1 1])
%     title('Proportion of Active Rich Neurons for Actual (Red) and Shuffled (Black) Data')
    subplot(2,4,ilen);
%     shadedErrorBar(AvalLen{ilen},Richness{ilen},Error{ilen},'r');hold on
    shadedErrorBar(AvalLen{ilen},rchnss{ilen},err{ilen},'r');hold on
    plot(AvalLen{ilen},.3*ones(1,length(AvalLen{ilen})),'k');
%     shadedErrorBar(AvalLen{ilen},rchnss_shuf_mean{ilen},err_shuf_mean{ilen},'k');
%     axis([1 len(ilen) .05 .8])
    axis([1 len(ilen) .1 .6])
    set(gca,'linewidth',2,'fontsize',16)
    xlabel(['avalanche length (x',num2str(BinSize),' ms)']);
% %     ylabel('mean # of active rich nodes')
end
% % save([dataset,'\AvalRichProfileShuf_AllDatasets_bs',num2str(BinSize),'ms.mat'],'AvalancheLen','RichnessShuff','ErrorShuff')
% % save([dataset,'\AvalRichProfile_AllDatasets_bs',num2str(BinSize),'ms.mat'],'AvalancheLen','Richness','Error')
% savefig([dataset,'\AvalRichActShuff_AllDatasets_bs',num2str(BinSize),'ms_subvivo'])
% saveas(H,[dataset,'\AvalRichActShuff_AllDatasets_bs',num2str(BinSize),'ms_subvivo.jpg'])


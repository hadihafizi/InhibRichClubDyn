clear all

dataSet = {'02-0','03-0','03-1','04-0','04-1','05-0','05-1','06-0',...
    '06-1','07-0','07-1','08-0','08-1','09-0','09-1'};

BinSize = 1; %ms

for i = 1:length(dataSet)
load([dataSet{i},'/asdf.mat'])

str_aval = Avalanche_detection(BinSize,asdf_raw);
thr = 45; % PDF significance level (threshold)
% load([dataset,'\PDF_1_16_30ms.mat']);
% sig_te = wgt.*PDF(:,:,thr); 
load([dataSet{i},'/PDF_NoSpur_thr45_',dataSet{i},'.mat'])
load([dataSet{i},'/wgts_1_16ms.mat']);
sig_te = wgt.*PDF; 
len = [2 5 8 12 15 20 25 30];
tic
[AvalLen,rchnss,err] = AvalRichProfileRev(len,str_aval,sig_te,0);
toc
%{
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
%}
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
    axis([1 len(ilen) .1 .5])
    set(gca,'linewidth',2,'fontsize',16)
    xlabel(['avalanche length (x',num2str(BinSize),' ms)']);
% %     ylabel('mean # of active rich nodes')
end

save([dataSet,'\AvalRichProfile_',dataSet{i},suffix,'_bs',num2str(BinSize),'ms.mat'],'AvalLen','rchnss','err')
savefig([dataSet,'\AvalRichActShuff_',dataSet{i},'_bs',num2str(BinSize),'ms',suffix])
saveas(H,[dataSet,'\AvalRichActShuff_',dataSet{i},'_bs',num2str(BinSize),'ms_',suffix,'.jpg'])
end
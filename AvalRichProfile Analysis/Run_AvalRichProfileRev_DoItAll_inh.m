clear all

dataSet = {'02-0','03-0','03-1','04-0','04-1','05-0','05-1','06-0',...
    '06-1','07-0','08-0','08-1','09-0','09-1'};
% dataSet = {'139', '151', '152', '163', '165', '168', '174'};
EIdatasets = [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25];

suffix = '_cor_inteRC';
% in-vitro
BinSize = 1; %ms
% in-vivo
% BinSize = 3; %ms
len = [2 5 8 12 15 20 25 30];
% parpool(2);
parfor i = 1:length(dataSet)

    PDF = load([dataSet{i},'/PDF_NoSpur_thr45_',dataSet{i},'.mat']);
    wgt = load([dataSet{i},'/wgts_1_16ms.mat']);
    sig_te = wgt.wgt.*PDF.PDF;
    
    exc_thresh = .9;
    inh_thresh = .9;
    % in-vitro
    neurontype = load(['C:\Users\mhafizi\Box Sync\Research\Beggs Lab\Projects\InhibID\Ian Stevenson lab\Neuron Type Analysis\Neuron Type Analysis\in vitro\DataSet', num2str(EIdatasets(i)), '.mat']);
    % in-vivo
%     neurontype = load(['C:\Users\Mahzad\Box Sync\Research\Beggs Lab\Projects\InhibID\Ian Stevenson lab\Neuron Type Analysis\Neuron Type Analysis\in vivo\unitactivity_Mouse ', num2str(dataSet{i}), '.mat']);
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

%     sig_te = sig_te(~ismember(1:length(sig_te), unidentified), ~ismember(1:length(sig_te), unidentified));
    oute = sum(sig_te,2);
    inte = sum(sig_te,1);
%     tote = oute + inte';
    [A1,B1] = sort(inte,'descend');
    num_rich = ceil(0.2*length(A1));
    r_idx = B1(1:num_rich);
    nr_idx = B1(num_rich+1:end);

    R_E = intersect(r_idx, ExcIdx, 'stable');
    R_I = intersect(r_idx, InhIdx, 'stable');
    NR_E = intersect(nr_idx, ExcIdx, 'stable');
    NR_I = intersect(nr_idx, InhIdx, 'stable');

%     ASDF = load([dataSet{i},'/asdf.mat']);
%     ASDF = load([dataSet{i},'/asdf_rest',dataSet{i},'_TimeCorOneMSGap.mat'])
%     asdf_raw = ASDF.asdf_raw(~ismember(1:length(ASDF.asdf_raw), unidentified));
%     asdf_raw{end}(1) = length(asdf_raw) - 2;
%     str_aval = Avalanche_detection(BinSize,asdf_raw);
%     parsave([dataSet{i},'/aval_list',dataSet{i},'_noUnID_bs',num2str(BinSize),'ms.mat'], {'str_aval'}, {str_aval});
    str_aval = load([dataSet{i},'/aval_list',dataSet{i},'_bs',num2str(BinSize),'ms.mat']);
    
%     temp = load([dataSet{i},'\AvalRichProfile_',dataSet{i},suffix,'_bs',num2str(BinSize),'ms.mat'])
%     AvalLen = temp.AvalLen;
%     rchnss = temp.rchnss;
%     err = temp.err;
    [AvalLen,rchnss,err] = AvalRichProfileRev_inh(len,str_aval.str_aval,r_idx);
    [~,rchnss_re,err_re] = AvalRichProfileRev_inh(len,str_aval.str_aval,R_E);
    [~,rchnss_ri,err_ri] = AvalRichProfileRev_inh(len,str_aval.str_aval,R_I);
    [~,rchnss_nr,err_nr] = AvalRichProfileRev_inh(len,str_aval.str_aval,nr_idx);
    [~,rchnss_nre,err_nre] = AvalRichProfileRev_inh(len,str_aval.str_aval,NR_E);
    [~,rchnss_nri,err_nri] = AvalRichProfileRev_inh(len,str_aval.str_aval,NR_I);

    %%
    %
    % plot the figures for all the chosen lengths
    H = figure(i);clf
    for ilen = 1 : length(len)
        set(H,'units','normalized','outerposition',[0 0 1 1])
        %     title('Proportion of Active Rich Neurons for Actual (Red) and Shuffled (Black) Data')
        subplot(2,4,ilen);
        %     shadedErrorBar(AvalLen{ilen},Richness{ilen},Error{ilen},'r');hold on
        shadedErrorBar(AvalLen{ilen},rchnss{ilen},err{ilen},'r');hold on
        shadedErrorBar(AvalLen{ilen},rchnss_re{ilen},err_re{ilen},'g');hold on
        shadedErrorBar(AvalLen{ilen},rchnss_ri{ilen},err_ri{ilen},'m');hold on
        shadedErrorBar(AvalLen{ilen},rchnss_nr{ilen},err_nr{ilen},'b');hold on
        shadedErrorBar(AvalLen{ilen},rchnss_nre{ilen},err_nre{ilen},'c');hold on
        shadedErrorBar(AvalLen{ilen},rchnss_nri{ilen},err_nri{ilen},'k');hold on
        plot(AvalLen{ilen},.3*ones(1,length(AvalLen{ilen})),'y');
        %     shadedErrorBar(AvalLen{ilen},rchnss_shuf_mean{ilen},err_shuf_mean{ilen},'k');
        %     axis([1 len(ilen) .05 .8])
        axis([1 len(ilen) 0 .7])
        set(gca,'linewidth',2,'fontsize',16)
        xlabel(['avalanche length (x',num2str(BinSize),' ms)']);
        % %     ylabel('mean # of active rich nodes')
    end
    
    parsave([dataSet{i},'\AvalRichProfile_',dataSet{i},suffix,'_bs',num2str(BinSize),'ms.mat'],{'AvalLen','rchnss','err'},{AvalLen,rchnss,err})
    parsave([dataSet{i},'\AvalRichProfile_',dataSet{i},suffix,'_bs',num2str(BinSize),'ms_RE.mat'],{'AvalLen','rchnss_re','err_re'},{AvalLen,rchnss_re,err_re})
    parsave([dataSet{i},'\AvalRichProfile_',dataSet{i},suffix,'_bs',num2str(BinSize),'ms_RI.mat'],{'AvalLen','rchnss_ri','err_ri'},{AvalLen,rchnss_ri,err_ri})
    parsave([dataSet{i},'\AvalRichProfile_',dataSet{i},suffix,'_bs',num2str(BinSize),'ms_NR.mat'],{'AvalLen','rchnss_nr','err_nr'},{AvalLen,rchnss_nr,err_nr})
    parsave([dataSet{i},'\AvalRichProfile_',dataSet{i},suffix,'_bs',num2str(BinSize),'ms_NRE.mat'],{'AvalLen','rchnss_nre','err_nre'},{AvalLen,rchnss_nre,err_nre})
    parsave([dataSet{i},'\AvalRichProfile_',dataSet{i},suffix,'_bs',num2str(BinSize),'ms_NRI.mat'],{'AvalLen','rchnss_nri','err_nri'},{AvalLen,rchnss_nri,err_nri})
    savefig(H,[dataSet{i},'\AvalRichActEI_',dataSet{i},'_bs',num2str(BinSize),'ms',suffix])
    saveas(H,[dataSet{i},'\AvalRichActEI_',dataSet{i},'_bs',num2str(BinSize),'ms_',suffix,'.jpg'])
end
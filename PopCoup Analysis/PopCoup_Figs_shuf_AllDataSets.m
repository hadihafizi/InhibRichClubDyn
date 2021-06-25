% Plot different distributions of similarity values for one data set at a
% time using the pruned network (no spurious connections).
%
% Required functions: None.
%
% Required files: PDF_1_16_30ms.mat
%                 PDF_NoSpur_thr[XX]_[dataSet]
%                 wgts_1_16ms.mat
%                 Sim_bin[X]ms.mat
%                 Sim_jit[X]bin_[X]ms.mat
%
% Hadi Hafizi, Apr. 2020
%
% Modified to use corrected (pruned for spurious connections) networks.
% HH, Aug. 2016
%
clear all

dataSets = {'02-0','03-0','03-1','04-0','04-1','05-0','05-1','06-0',...
    '06-1','07-0','07-1','08-0','08-1','09-0','09-1'};
% dataSets = {'139', '151', '152', '163', '165', '168', '174'};

thr = 45;
bs = 1; % bin size [ms]
JitBin = 30;
% load asdf.mat;

%%

Coup_r_r_norm = [];
Coup_r_nr_norm = [];
Coup_nr_r_norm = [];
Coup_nr_nr_norm = [];
Coup_r_r_shuf_norm = [];
Coup_r_nr_shuf_norm = [];
Coup_nr_r_shuf_norm = [];
Coup_nr_nr_shuf_norm = [];

for i = 1:length(dataSets)
    load(['PopCoupAll_',dataSets{i},'_cor_toteRC.mat'])
%     load(['PopCoupAll_',dataSets{i},'_cor_toteRC_bs',num2str(bs),'ms.mat'])

    [A1, B1] = size(Coup_r_r);[A2, B2] = size(Coup_nr_nr);
    interval = (B1-1)/2;
    TLag = -interval:interval;
    interval_shuf = 50;
    % normalize for the number of nodes and subtract baseline.
    for ii = 1:A1
        Coup_r_r(ii,:) = (Coup_r_r(ii,:)-Coup_r_r(ii,1))./(A1-1);
        Coup_r_nr(ii,:) = (Coup_r_nr(ii,:)-Coup_r_nr(ii,1))./(A2-1);
%         Coup_r_r(ii,:) = (Coup_r_r(ii,:))./(A1-1);
%         Coup_r_nr(ii,:) = (Coup_r_nr(ii,:))./(A2-1);
    end
    for ii = 1:A2
        Coup_nr_nr(ii,:) = (Coup_nr_nr(ii,:)-Coup_nr_nr(ii,1))./(A2-1);
        Coup_nr_r(ii,:) = (Coup_nr_r(ii,:)-Coup_nr_r(ii,1))./(A1-1);
%         Coup_nr_nr(ii,:) = (Coup_nr_nr(ii,:))./(A2-1);
%         Coup_nr_r(ii,:) = (Coup_nr_r(ii,:))./(A1-1);
    end
    
    Coup_r_r_norm = [Coup_r_r_norm;Coup_r_r];
    Coup_r_nr_norm = [Coup_r_nr_norm;Coup_r_nr];
    Coup_nr_r_norm = [Coup_nr_r_norm;Coup_nr_r];
    Coup_nr_nr_norm = [Coup_nr_nr_norm;Coup_nr_nr];
    
    Coup_r_r_shuf_sum = zeros(A1,2*interval_shuf+1);
    Coup_nr_nr_shuf_sum = zeros(A2,2*interval_shuf+1);
    Coup_r_nr_shuf_sum = zeros(A1,2*interval_shuf+1);
    Coup_nr_r_shuf_sum = zeros(A2,2*interval_shuf+1);
    
    for ii=1:20
        load(['PopCoupAll_', dataSets{i}, '_cor_toteRC_bs',num2str(bs),'ms_shuf', num2str(ii), '.mat'])
        
        Coup_r_r_shuf_sum = Coup_r_r + Coup_r_r_shuf_sum;
        Coup_nr_nr_shuf_sum = Coup_nr_nr + Coup_nr_nr_shuf_sum;
        Coup_r_nr_shuf_sum = Coup_r_nr + Coup_r_nr_shuf_sum;
        Coup_nr_r_shuf_sum = Coup_nr_r + Coup_nr_r_shuf_sum;
        
    end
    Coup_r_r_shuf_sum = Coup_r_r_shuf_sum/ii;
    Coup_nr_nr_shuf_sum = Coup_nr_nr_shuf_sum/ii;
    Coup_r_nr_shuf_sum = Coup_r_nr_shuf_sum/ii;
    Coup_nr_r_shuf_sum = Coup_nr_r_shuf_sum/ii;

    % normalize for the number of nodes and subtract baseline.
   for ii = 1:A1
        Coup_r_r_shuf(ii,:) = (Coup_r_r_shuf_sum(ii,:)-Coup_r_r_shuf_sum(ii,1))./(A1-1);
        Coup_r_nr_shuf(ii,:) = (Coup_r_nr_shuf_sum(ii,:)-Coup_r_nr_shuf_sum(ii,1))./(A2-1);
%         Coup_r_r_shuf(ii,:) = (Coup_r_r_shuf_sum(ii,:))./(A1-1);
%         Coup_r_nr_shuf(ii,:) = (Coup_r_nr_shuf_sum(ii,:))./(A2-1);
    end
    for ii = 1:A2
        Coup_nr_nr_shuf(ii,:) = (Coup_nr_nr_shuf_sum(ii,:)-Coup_nr_nr_shuf_sum(ii,1))./(A2-1);
        Coup_nr_r_shuf(ii,:) = (Coup_nr_r_shuf_sum(ii,:)-Coup_nr_r_shuf_sum(ii,1))./(A1-1);
%         Coup_nr_nr_shuf(ii,:) = (Coup_nr_nr_shuf_sum(ii,:))./(A2-1);
%         Coup_nr_r_shuf(ii,:) = (Coup_nr_r_shuf_sum(ii,:))./(A1-1);
    end
    
    Coup_r_r_shuf_norm = [Coup_r_r_shuf_norm;Coup_r_r_shuf];
    Coup_r_nr_shuf_norm = [Coup_r_nr_shuf_norm;Coup_r_nr_shuf];
    Coup_nr_r_shuf_norm = [Coup_nr_r_shuf_norm;Coup_nr_r_shuf];
    Coup_nr_nr_shuf_norm = [Coup_nr_nr_shuf_norm;Coup_nr_nr_shuf];
    
end

PopCoupAll = [Coup_r_r_norm(:,interval+1); Coup_nr_r_norm(:,interval+1); ...
    Coup_nr_nr_norm(:,interval+1); Coup_r_nr_norm(:,interval+1); ...
    Coup_r_r_shuf_norm(:,interval_shuf+1); Coup_nr_r_shuf_norm(:,interval_shuf+1); ...
    Coup_nr_nr_shuf_norm(:,interval_shuf+1); Coup_r_nr_shuf_norm(:,interval_shuf+1)];

[~, edges] = histcounts(PopCoupAll);
[rrcount, ~] = histcounts(Coup_r_r_norm(:,interval+1),edges);
[nrrcount, ~] = histcounts(Coup_nr_r_norm(:,interval+1),edges);
[nrnrcount, ~] = histcounts(Coup_nr_nr_norm(:,interval+1),edges);
[rnrcount, ~] = histcounts(Coup_r_nr_norm(:,interval+1),edges);
[rrcount_shuf, ~] = histcounts(Coup_r_r_shuf_norm(:,interval_shuf+1),edges);
[nrrcount_shuf, ~] = histcounts(Coup_nr_r_shuf_norm(:,interval_shuf+1),edges);
[nrnrcount_shuf, ~] = histcounts(Coup_nr_nr_shuf_norm(:,interval_shuf+1),edges);
[rnrcount_shuf, ~] = histcounts(Coup_r_nr_shuf_norm(:,interval_shuf+1),edges);

rrcount = rrcount/sum(rrcount);
nrnrcount = nrnrcount/sum(nrnrcount);
rnrcount = rnrcount/sum(rnrcount);
nrrcount = nrrcount/sum(nrrcount);
rrcount_shuf = rrcount_shuf/sum(rrcount_shuf);
nrnrcount_shuf = nrnrcount_shuf/sum(nrnrcount_shuf);
rnrcount_shuf = rnrcount_shuf/sum(rnrcount_shuf);
nrrcount_shuf = nrrcount_shuf/sum(nrrcount_shuf);

centers = movsum(edges, 2)/2;
centers = centers(2:end);

%% t-test

ttest(Coup_r_r_norm(:,interval+1), Coup_r_r_shuf_norm(:,interval_shuf+1))



%%

figure;
subplot(2,2,1);
plot(centers,rrcount,'r--.','MarkerSize',10); hold on
plot(centers,rrcount_shuf,'y--.','MarkerSize',10);
legend('RR','RR_{shuf}');
axis([-100 500 0 .9]);
title('R-->R', 'FontSize',16)
xlabel('Population Coupling [au]');
ylabel('Probability');
set(gca,'tickDir','out','color','none','fontname','arial','linewidth',1,'fontsize',12); box off;

subplot(2,2,2);
plot(centers,nrnrcount,'k--.','MarkerSize',10); hold on
plot(centers,nrnrcount_shuf,'y--.','MarkerSize',10);
legend('NRNR','NRNR_{shuf}');
axis([-100 500 0 .9]);
title('NR-->NR', 'FontSize',16)
xlabel('Population Coupling [au]');
ylabel('Probability');
set(gca,'tickDir','out','color','none','fontname','arial','linewidth',1,'fontsize',12); box off;

subplot(2,2,3);
plot(centers,rnrcount,'b--.','MarkerSize',10); hold on
plot(centers,rnrcount_shuf,'y--.','MarkerSize',10);
legend('RNR','RNR_{shuf}');
axis([-100 500 0 .9]);
title('R-->NR', 'FontSize',16)
xlabel('Population Coupling [au]');
ylabel('Probability');
set(gca,'tickDir','out','color','none','fontname','arial','linewidth',1,'fontsize',12); box off;

subplot(2,2,4);
plot(centers,nrrcount,'c--.','MarkerSize',10); hold on
plot(centers,nrrcount_shuf,'y--.','MarkerSize',10);
legend('NRR','NRR_{shuf}');
axis([-100 500 0 .9]);
title('NR-->R', 'FontSize',16)
xlabel('Population Coupling [au]');
ylabel('Probability');
set(gca,'tickDir','out','color','none','fontname','arial','linewidth',1,'fontsize',12); box off;

% legend('R-R','NR-NR','R-NR','NR-R','R-R_shuf','NR-NR_shuf','R-NR_shuf','NR-R_shuf');

%% KL-Div PopCoup

KL_RRvsNRNR = zeros(length(dataSets),1); KL_RRvsRNR = zeros(length(dataSets),1); 
KL_NRNRvsRNR = zeros(length(dataSets),1); KL_RRvsNRR = zeros(length(dataSets),1); 
KL_NRNRvsNRR = zeros(length(dataSets),1); KL_RNRvsNRR = zeros(length(dataSets),1); 
for i = 1:length(dataSets)
    load(['PopCoupAll_',dataSets{i},'_cor_toteRC.mat'])
    
    [A1, B1] = size(Coup_r_r);[A2, B2] = size(Coup_nr_nr);
    interval = (B1 - 1)/2;
    
    for ii = 1:A1
        Coup_r_r_norm(ii,:) = (Coup_r_r(ii,:)-Coup_r_r(ii,1))./(A1-1);
        Coup_r_nr_norm(ii,:) = (Coup_r_nr(ii,:)-Coup_r_nr(ii,1))./(A2-1);
    end
    for ii = 1:A2
        Coup_nr_nr_norm(ii,:) = (Coup_nr_nr(ii,:)-Coup_nr_nr(ii,1))./(A2-1);
        Coup_nr_r_norm(ii,:) = (Coup_nr_r(ii,:)-Coup_nr_r(ii,1))./(A1-1);
    end
    
    PopCoupAll = [Coup_r_r_norm(:,interval+1); Coup_nr_r_norm(:,interval+1); ...
        Coup_nr_nr_norm(:,interval+1); Coup_r_nr_norm(:,interval+1)];
    [~, edges] = histcounts(PopCoupAll);
    [rrcount, ~] = histcounts(Coup_r_r_norm(:,interval+1),edges);
    [nrrcount, ~] = histcounts(Coup_nr_r_norm(:,interval+1),edges);
    [nrnrcount, ~] = histcounts(Coup_nr_nr_norm(:,interval+1),edges);
    [rnrcount, ~] = histcounts(Coup_r_nr_norm(:,interval+1),edges);
    
    rrcount = rrcount/sum(rrcount);
    nrnrcount = nrnrcount/sum(nrnrcount);
    rnrcount = rnrcount/sum(rnrcount);
    nrrcount = nrrcount/sum(nrrcount);
    
    centers = movsum(edges, 2)/2;
    centers = centers(2:end);
    
    KL_RRvsNRNR(i) = kldiv(centers,rrcount/sum(rrcount) + eps, nrnrcount/sum(nrnrcount) + eps);
    KL_RRvsRNR(i) = kldiv(centers,rrcount/sum(rrcount) + eps, rnrcount/sum(rnrcount) + eps);
    KL_RRvsNRR(i) = kldiv(centers,rrcount/sum(rrcount) + eps, nrrcount/sum(nrrcount) + eps);
    KL_NRNRvsRNR(i) = kldiv(centers,rnrcount/sum(rnrcount) + eps, nrnrcount/sum(nrnrcount) + eps);
    KL_NRNRvsNRR(i) = kldiv(centers,nrrcount/sum(nrrcount) + eps, nrnrcount/sum(nrnrcount) + eps);
    KL_RNRvsNRR(i) = kldiv(centers,nrrcount/sum(nrrcount) + eps, rnrcount/sum(rnrcount) + eps);

end

kldiv_all = [KL_RRvsNRNR, KL_RRvsRNR, KL_RRvsNRR, KL_NRNRvsRNR, KL_NRNRvsNRR, KL_RNRvsNRR];
figure;
bar([1:3, 5:7], mean(kldiv_all)); hold on
xticklabels({'D_{KL}(RR||NRNR)','D_{KL}(RR||RNR)','D_{KL}(RR||NRR)','D_{KL}(NRNR||RNR)','D_{KL}(RNR||NRR)','D_{KL}(RNR||NRR)'});
ylabel('KL-Divergence [bits]')
er = errorbar([1:3, 5:7], mean(kldiv_all), std(kldiv_all)/sqrt(15));
er.Color = [0 0 0];
er.LineStyle = 'none';
set(gca,'tickDir','out','color','none','fontname','arial','linewidth',1,'fontsize',9); box off;

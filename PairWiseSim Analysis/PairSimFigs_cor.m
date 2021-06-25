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
% Hadi Hafizi, Nov. 2015
%
% Modified to use corrected (pruned for spurious connections) networks.
% HH, Aug. 2016
%
clear all

dataSet = '020';
thr = 45;
JitBin = 30;
% load asdf.mat;  
load([dataSet,'\PDF_1_16_30ms.mat']); 
load([dataSet,'\wgts_1_16ms.mat'])
W = PDF(:,:,45).*wgt;oute = sum(W,2);
[Rich Rich_idx] = sort(oute,'descend');
num_rich = ceil(0.2*length(Rich));
Rich_idx = Rich_idx(1:num_rich); nonRich_idx = setdiff(1:length(Rich),Rich_idx);
nr_oute = oute(nonRich_idx);
[nr_vals, nr_idx] = sort(nr_oute,'descend');

load([dataSet,'\Sim_bin1ms.mat'])
PairSim = S;
clear S
S_r = PairSim(1:length(Rich_idx),1:length(Rich_idx));
S_nr = PairSim(length(Rich_idx)+1:end,length(Rich_idx)+1:end);
S_r_nr = PairSim(1:length(Rich_idx),length(Rich_idx)+1:end);


S_shuf_mean = mean(S_shuf,3);
S_shuf_r = mean(S_shuf(1:length(Rich_idx),1:length(Rich_idx),:),3);
S_shuf_nr = mean(S_shuf(length(Rich_idx)+1:end,length(Rich_idx)+1:end,:),3);
S_shuf_r_nr = mean(S_shuf(1:length(Rich_idx),length(Rich_idx)+1:end,:),3);

load([dataSet,'\Sim_jit',num2str(JitBin),'bin_1ms.mat'])
SimJit = S;
clear S

SimJit = mean(SimJit,3);
SimJit_r = SimJit(1:length(Rich_idx),1:length(Rich_idx));
SimJit_nr = SimJit(length(Rich_idx)+1:end,length(Rich_idx)+1:end);
SimJit_r_nr = SimJit(1:length(Rich_idx),length(Rich_idx)+1:end);

% re-arrange for new RC members
[~, Rich_idx] = sort(oute,'descend');
[~, new_sort] = sort(Rich_idx,'ascend');
PairSim_cor = PairSim(new_sort,new_sort);
S_shuf_mean_cor = S_shuf_mean(new_sort,new_sort);
SimJit_cor = SimJit(new_sort,new_sort);
load([dataSet,'\PDF_NoSpur_thr',num2str(thr),'_',dataSet])
W = PDF_cor.*wgt;oute = sum(W,2);
[Rich, Rich_idx_cor] = sort(oute,'descend');
num_rich = ceil(0.2*length(Rich));
nonRich_idx_cor = Rich_idx_cor(num_rich+1:end);
Rich_idx_cor = Rich_idx_cor(1:num_rich);

S_r = PairSim_cor(Rich_idx_cor,Rich_idx_cor);
S_nr = PairSim_cor(nonRich_idx_cor,nonRich_idx_cor);
S_r_nr = PairSim_cor(Rich_idx_cor,nonRich_idx_cor);

S_shuf_r = S_shuf_mean_cor(Rich_idx_cor,Rich_idx_cor);
S_shuf_nr = S_shuf_mean_cor(nonRich_idx_cor,nonRich_idx_cor);
S_shuf_r_nr = S_shuf_mean_cor(Rich_idx_cor,nonRich_idx_cor);

SimJit_r = SimJit(Rich_idx_cor,Rich_idx_cor);
SimJit_nr = SimJit(nonRich_idx_cor,nonRich_idx_cor);
SimJit_r_nr = SimJit(Rich_idx_cor,nonRich_idx_cor);

%%
nbin = 50;
maxSim = max([max(S_r(:,:)),max(S_nr(:,:)),...
    max(S_r_nr(:,:))]);
minSim = min([min(S_r(:,:)),min(S_nr(:,:)),...
    min(S_r_nr(:,:))]);
bnedg = minSim:((maxSim-minSim)/nbin):maxSim;

[rrcount, rrcenter] = hist(reshape(S_r(triu(ones(length(Rich_idx_cor)),1)...
    ==ones(length(Rich_idx_cor))),1,length(Rich_idx_cor)*(length(Rich_idx_cor)-1)/2),...
    bnedg); 
[nrnrcount, nrnrcenter] = hist(reshape(S_nr(triu(ones(length(nonRich_idx_cor)),1)...
    ==ones(length(nonRich_idx_cor))),1,length(nonRich_idx_cor)*(length(nonRich_idx_cor)-1)/2),...
    bnedg); 
[nrrcount, nrrcenter] = hist(reshape(S_r_nr,1,size(S_r_nr,1)*size(S_r_nr,2)),bnedg); 

rrcount = rrcount/sum(rrcount);
nrnrcount = nrnrcount/sum(nrnrcount);
nrrcount = nrrcount/sum(nrrcount);

figure(107);
%{
frr = fit(rrcenter',rrcount','gauss1');
fnrnr = fit(nrnrcenter',nrnrcount','gauss1');
fnrr = fit(nrrcenter',nrrcount','gauss1');

plot(frr,rrcenter,rrcount); hold on
plot(fnrnr,nrnrcenter,nrnrcount); hold on
plot(fnrr,nrrcenter,nrrcount); hold on
%}
plot(rrcenter,rrcount,'r--.','MarkerSize',10); hold on
plot(nrnrcenter,nrnrcount,'b--.','MarkerSize',10); hold on
plot(nrrcenter,nrrcount,'k--.','MarkerSize',10); hold on

% legend('','R-R','','NR-NR','','R-NR');
legend('R-R','NR-NR','R-NR');
title('Jaccard Coefficients Distributions', 'FontSize',16)
xlabel('Jaccard Coefficient');
ylabel('Probability');
set(gca, 'FontSize',16,'xscale','log');

%%
%{
nbin = 50;
maxSim = max([max(S_r(:,:)),max(S_shuf_r(:,:)),...
    max(SimJit_r(:,:))]);
minSim = min([min(S_r(:,:)),min(S_shuf_r(:,:)),...
    min(SimJit_r(:,:))]);
bnedg = minSim:((maxSim-minSim)/nbin):maxSim;
%}
% S_r % S_nr % S_r_nr % S_shuf_r % S_shuf_nr % S_shuf_r_nr % SimJit_r
% SimJit_nr % SimJit_r_nr
[rrcount, rrcenter] = hist(reshape(S_r(triu(ones(length(Rich_idx_cor)),1)...
    ==ones(length(Rich_idx_cor))),1,length(Rich_idx_cor)*(length(Rich_idx_cor)-1)/2),...
    bnedg); 
[nrnrcount, nrnrcenter] = hist(reshape(S_shuf_r(triu(ones(length(Rich_idx_cor)),1)...
    ==ones(length(Rich_idx_cor))),1,length(Rich_idx_cor)*(length(Rich_idx_cor)-1)/2),...
    bnedg); 
[nrrcount, nrrcenter] = hist(reshape(SimJit_r(triu(ones(length(Rich_idx_cor)),1)...
    ==ones(length(Rich_idx_cor))),1,length(Rich_idx_cor)*(length(Rich_idx_cor)-1)/2),...
    bnedg); 

rrcount = rrcount/sum(rrcount);
nrnrcount = nrnrcount/sum(nrnrcount);
nrrcount = nrrcount/sum(nrrcount);

figure(108);
%{
frr = fit(rrcenter',rrcount','gauss1');
fnrnr = fit(nrnrcenter',nrnrcount','gauss1');
fnrr = fit(nrrcenter',nrrcount','gauss1');

plot(frr,rrcenter,rrcount); hold on
plot(fnrnr,nrnrcenter,nrnrcount); hold on
plot(fnrr,nrrcenter,nrrcount); hold on
%}
plot(rrcenter,rrcount,'r--.','MarkerSize',10); hold on
plot(nrnrcenter,nrnrcount,'b--.','MarkerSize',10); hold on
plot(nrrcenter,nrrcount,'k--.','MarkerSize',10); hold on

% legend('','Actual','','Shuffled','','Jittered');
legend('Actual','Shuffled','Jittered');
title('Jaccard Coefficients Distributions', 'FontSize',16)
xlabel('Jaccard Coefficient');
ylabel('Probability');
set(gca, 'FontSize',16,'xscale','log');

%%

rrcount_rand = zeros(size(rrcount));
nrnrcount_rand = zeros(size(nrnrcount));
nrrcount_rand = zeros(size(nrrcount));
for i = 1:1000
    randIdx20 = randperm(size(oute,1),ceil(.2*size(oute,1)));
    randIdx80 = setdiff(1:size(oute,1),randIdx20);
    S_r_rand = PairSim_cor(randIdx20,randIdx20);
    [rrcount, rrcenter] = hist(reshape(S_r_rand(triu(ones(length(Rich_idx_cor)),1)...
        ==ones(length(Rich_idx_cor))),1,length(Rich_idx_cor)*(length(Rich_idx_cor)-1)/2),...
        bnedg);
    rrcount_rand = rrcount + rrcount_rand;
    S_nr_rand = PairSim_cor(randIdx80,randIdx80);
    [nrnrcount, nrnrcenter] = hist(reshape(S_nr_rand(triu(ones(length(nonRich_idx_cor)),1)...
        ==ones(length(nonRich_idx_cor))),1,length(nonRich_idx_cor)*(length(nonRich_idx_cor)-1)/2),...
        bnedg);
    nrnrcount_rand = nrnrcount + nrnrcount_rand;
    S_r_nr_rand = PairSim_cor(randIdx20,randIdx80);
    [nrrcount, nrrcenter] = hist(reshape(S_r_nr_rand,1,size(S_r_nr,1)*size(S_r_nr,2)),bnedg);
    nrrcount_rand = nrrcount + nrrcount_rand;
end

rrcount_rand = rrcount_rand/sum(rrcount_rand);
nrnrcount_rand = nrnrcount_rand/sum(nrnrcount_rand);
nrrcount_rand = nrrcount_rand/sum(nrrcount_rand);

figure(107); hold on
plot(rrcenter,rrcount_rand,'g--.','MarkerSize',10); hold on
plot(rrcenter,nrnrcount_rand,'y--.','MarkerSize',10); hold on
plot(rrcenter,nrrcount_rand,'c--.','MarkerSize',10); hold on
legend('R-R','NR-NR','R-NR','rnd 20%','rnd 80%','btw rnd');


%%

F1 = figure;
subplot(1,3,1);imagesc(PairSim_cor([Rich_idx_cor;nonRich_idx_cor],[Rich_idx_cor;nonRich_idx_cor])+...
    PairSim_cor([Rich_idx_cor;nonRich_idx_cor],[Rich_idx_cor;nonRich_idx_cor])');colorbar;colormap jet;
title('PairSimAct');
set(gca,'FontSize',16);
subplot(1,3,2);imagesc(SimJit_cor([Rich_idx_cor;nonRich_idx_cor],[Rich_idx_cor;nonRich_idx_cor])+...
    SimJit_cor([Rich_idx_cor;nonRich_idx_cor],[Rich_idx_cor;nonRich_idx_cor])');colorbar;colormap jet;
title(['PairSimJit ',num2str(JitBin),'ms']);
set(gca,'FontSize',16);
subplot(1,3,3);imagesc(S_shuf_mean_cor([Rich_idx_cor;nonRich_idx_cor],[Rich_idx_cor;nonRich_idx_cor])+...
    S_shuf_mean_cor([Rich_idx_cor;nonRich_idx_cor],[Rich_idx_cor;nonRich_idx_cor])');colorbar;colormap jet;
title('PairSimShuff');
set(gca,'FontSize',16);
% cd Figures
% savefig(['PairSimMatActShuffJit',num2str(JitBin),'ms'])
% saveas(F1,['PairSimMatActShuffJit',num2str(JitBin),'ms.jpg'])
% cd ..

F2 = figure;
NumBins = 90;
h1 = histogram(triu(S_r,1),NumBins,'normalization','probability');
bnedg = h1.BinEdges;hold on
histogram(S_r_nr,bnedg,'normalization','probability')
hold on
histogram(triu(S_nr,1),bnedg,'normalization','probability')
legend('Rich','r-nr','nonRich')
title('Pairwise Similarity Histogram','FontSize',16);
xlabel('Similarity','FontSize',16);
ylabel('Probability','FontSize',16);
set(gca,'FontSize',16);
axis([0 0.035 -inf inf])
% cd Figures
% savefig('PairSimHistAct_all')
% saveas(F2,'PairSimHistAct_all.jpg')
% cd ..

F3 = figure;
% NumBins = 50;
h2 = histogram(triu(S_r,1),NumBins,'normalization','probability');
bnedg = h2.BinEdges;hold on
histogram(triu(SimJit_r,1),bnedg,'normalization','probability')
hold on
histogram(triu(S_shuf_r,1),bnedg,'normalization','probability')
legend('ActualRich',['JitRich ',num2str(JitBin),'ms'],'ShuffRich')
title('Pairwise Similarity Histogram','FontSize',16);
xlabel('Similarity','FontSize',16);
ylabel('Probability','FontSize',16);
set(gca,'FontSize',16);
axis([0 0.035 -inf inf])
% cd Figures
% savefig(['PairSimHistActShuffJit',num2str(JitBin),'ms_r'])
% saveas(F3,['PairSimHistActShuffJit',num2str(JitBin),'ms_r.jpg'])
% cd ..

F4 = figure;
% NumBins = 50;
h3 = histogram(S_r_nr,NumBins,'normalization','probability');
bnedg = h3.BinEdges;hold on
histogram(SimJit_r_nr,bnedg,'normalization','probability')
hold on
histogram(S_shuf_r_nr,bnedg,'normalization','probability')
legend('Actual-r-nr',['Jit-r-nr ',num2str(JitBin),'ms'],'Shuff-r-nr')
title('Pairwise Similarity Histogram','FontSize',16);
xlabel('Similarity','FontSize',16);
ylabel('Probability','FontSize',16);
set(gca,'FontSize',16);
axis([0 0.035 -inf inf])
% cd Figures
% savefig(['PairSimHistActShuffJit',num2str(JitBin),'ms_r_nr'])
% saveas(F4,['PairSimHistActShuffJit',num2str(JitBin),'ms_r_nr.jpg'])
% cd ..

F5 = figure;
% NumBins = 50;
h4 = histogram(triu(PairSim_cor,1),NumBins,'normalization','probability');
bnedg = h4.BinEdges;hold on
histogram(triu(SimJit_cor,1),bnedg,'normalization','probability')
hold on
histogram(triu(S_shuf_mean_cor,1),bnedg,'normalization','probability')
legend('Actual',['Jitter ',num2str(JitBin),'ms'],'Shuff')
title('Pairwise Similarity Histogram','FontSize',16);
xlabel('Similarity','FontSize',16);
ylabel('Probability','FontSize',16);
set(gca,'FontSize',16);
axis([0 0.035 -inf inf])
% cd Figures
% savefig(['PairSimHistActShuffJit',num2str(JitBin),'ms_all'])
% saveas(F4,['PairSimHistActShuffJit',num2str(JitBin),'ms_all.jpg'])
% cd ..

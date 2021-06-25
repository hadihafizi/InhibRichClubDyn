
clear all

dataSets = {'02-0','03-0','03-1', '04-0', '04-1','05-0', '05-1','06-0',...
    '06-1','07-0','07-1','08-0', '08-1','09-0','09-1'};
Coup_r_r_norm = [];
Coup_r_nr_norm = [];
Coup_nr_r_norm = [];
Coup_nr_nr_norm = [];
Coup_r_r_shuf_norm = [];
Coup_r_nr_shuf_norm = [];
Coup_nr_r_shuf_norm = [];
Coup_nr_nr_shuf_norm = [];
for i = 1:length(dataSets)
    load(['PopCoupAll_', dataSets{i}, '_cor_toteRC.mat'])
    [Coup_r_r,Coup_r_nr,Coup_nr_r,Coup_nr_nr,...
        Coup_r_r_shuf,Coup_r_nr_shuf,Coup_nr_r_shuf,Coup_nr_nr_shuf] = PopCoup_Figs(Coup_r_r, Coup_r_nr, Coup_nr_r, Coup_nr_nr, dataSets{i}, 0);
    Coup_r_r_norm = [Coup_r_r_norm;Coup_r_r];
    Coup_r_nr_norm = [Coup_r_nr_norm;Coup_r_nr];
    Coup_nr_r_norm = [Coup_nr_r_norm;Coup_nr_r];
    Coup_nr_nr_norm = [Coup_nr_nr_norm;Coup_nr_nr];
    Coup_r_r_shuf_norm = [Coup_r_r_shuf_norm;Coup_r_r_shuf];
    Coup_r_nr_shuf_norm = [Coup_r_nr_shuf_norm;Coup_r_nr_shuf];
    Coup_nr_r_shuf_norm = [Coup_nr_r_shuf_norm;Coup_nr_r_shuf];
    Coup_nr_nr_shuf_norm = [Coup_nr_nr_shuf_norm;Coup_nr_nr_shuf];
end

%% Actual Data
nbin = 100;
maxCoup = max([max(Coup_r_r_norm(:,:)),max(Coup_nr_r_norm(:,:)),...
    max(Coup_nr_nr_norm(:,:)),max(Coup_r_nr_norm(:,:))]);
minCoup = min([min(Coup_r_r_norm(:,:)),min(Coup_nr_r_norm(:,:)),...
    min(Coup_nr_nr_norm(:,:)),min(Coup_r_nr_norm(:,:))]);
bnedg = minCoup:((maxCoup-minCoup)/nbin):maxCoup;
%%
[rrcount, rrcenter] = hist(Coup_r_r_norm(:,51),bnedg); 
[nrrcount, nrrcenter] = hist(Coup_nr_r_norm(:,51),bnedg); 
[nrnrcount, nrnrcenter] = hist(Coup_nr_nr_norm(:,51),bnedg); 
[rnrcount, rnrcenter] = hist(Coup_r_nr_norm(:,51),bnedg); 

rrcount = rrcount/sum(rrcount);
nrnrcount = nrnrcount/sum(nrnrcount);
rnrcount = rnrcount/sum(rnrcount);
nrrcount = nrrcount/sum(nrrcount);

figure;
plot(rrcenter,rrcount,'r--.','MarkerSize',10); hold on
plot(nrnrcenter,nrnrcount,'k--.','MarkerSize',10); hold on
plot(rnrcenter,rnrcount,'b--.','MarkerSize',10); hold on
plot(nrrcenter,nrrcount,'c--.','MarkerSize',10); hold on

legend('R-R','NR-NR','R-NR','NR-R');
title('Population Coupling Distributions', 'FontSize',16)
xlabel('Population Coupling [au]');
ylabel('Probability');
set(gca, 'FontSize',16);

%%
ax = figure(); hold on
rr = histogram(Coup_r_r_norm(:,51),bnedg, 'normalization', 'probability', 'DisplayStyle','stairs', 'LineWidth', 2);
nrnr = histogram(Coup_nr_nr_norm(:,51),bnedg, 'normalization', 'probability', 'DisplayStyle','stairs', 'LineWidth', 2);
rnr = histogram(Coup_r_nr_norm(:,51),bnedg, 'normalization', 'probability', 'DisplayStyle','stairs', 'LineWidth', 2);
nrr = histogram(Coup_nr_r_norm(:,51),bnedg, 'normalization', 'probability', 'DisplayStyle','stairs', 'LineWidth', 2);

legend('R-R','NR-NR','R-NR','NR-R');
title('Population Coupling Distributions', 'FontSize',16)
xlabel('Population Coupling [au]');
ylabel('Probability Density');
set(gca, 'FontSize',16);

figure;
plot(rr.BinEdges(2:end), rr.Values,'r--.','MarkerSize',10); hold on
plot(nrnr.BinEdges(2:end), nrnr.Values,'k--.','MarkerSize',10); hold on
plot(rnr.BinEdges(2:end), rnr.Values,'b--.','MarkerSize',10); hold on
plot(nrr.BinEdges(2:end), nrr.Values,'c--.','MarkerSize',10);

legend('R-R','NR-NR','R-NR','NR-R');
title('Population Coupling Distributions', 'FontSize',16)
xlabel('Population Coupling [au]');
ylabel('Probability Density');
set(gca, 'FontSize',16);

%% Shuffled Data
nbin = 100;
maxCoup = max([max(Coup_r_r_shuf_norm(:,:)),max(Coup_nr_r_shuf_norm(:,:)),...
    max(Coup_nr_nr_shuf_norm(:,:)),max(Coup_r_nr_shuf_norm(:,:))]);
minCoup = min([min(Coup_r_r_shuf_norm(:,:)),min(Coup_nr_r_shuf_norm(:,:)),...
    min(Coup_nr_nr_shuf_norm(:,:)),min(Coup_r_nr_shuf_norm(:,:))]);
bnedg = minCoup:((maxCoup-minCoup)/nbin):maxCoup;
%%
[rrcount, rrcenter] = hist(Coup_r_r_shuf_norm(:,51),bnedg); 
[nrrcount, nrrcenter] = hist(Coup_nr_r_shuf_norm(:,51),bnedg); 
[nrnrcount, nrnrcenter] = hist(Coup_nr_nr_shuf_norm(:,51),bnedg); 
[rnrcount, rnrcenter] = hist(Coup_r_nr_shuf_norm(:,51),bnedg); 

rrcount = rrcount/sum(rrcount);
nrnrcount = nrnrcount/sum(nrnrcount);
rnrcount = rnrcount/sum(rnrcount);
nrrcount = nrrcount/sum(nrrcount);

figure;
plot(rrcenter,rrcount,'r--.','MarkerSize',10); hold on
plot(nrnrcenter,nrnrcount,'k--.','MarkerSize',10); hold on
plot(rnrcenter,rnrcount,'b--.','MarkerSize',10); hold on
plot(nrrcenter,nrrcount,'c--.','MarkerSize',10); hold on

legend('R-R','NR-NR','R-NR','NR-R');
title('Population Coupling Distributions-Shuffled', 'FontSize',16)
xlabel('Population Coupling [au]');
ylabel('Probability');
set(gca, 'FontSize',16);

%%
%%
ax = figure(); hold on
rr = histogram(Coup_r_r_shuf_norm(:,51),bnedg, 'normalization', 'probability', 'DisplayStyle','stairs', 'LineWidth', 2);
nrnr = histogram(Coup_nr_nr_shuf_norm(:,51),bnedg, 'normalization', 'probability', 'DisplayStyle','stairs', 'LineWidth', 2);
rnr = histogram(Coup_r_nr_shuf_norm(:,51),bnedg, 'normalization', 'probability', 'DisplayStyle','stairs', 'LineWidth', 2);
nrr = histogram(Coup_nr_r_shuf_norm(:,51),bnedg, 'normalization', 'probability', 'DisplayStyle','stairs', 'LineWidth', 2);

legend('R-R','NR-NR','R-NR','NR-R');
title('Population Coupling Distributions', 'FontSize',16)
xlabel('Population Coupling [au]');
ylabel('Probability Density');
set(gca, 'FontSize',16);

figure;
plot(rr.BinEdges(2:end), rr.Values,'r--.','MarkerSize',10); hold on
plot(nrnr.BinEdges(2:end), nrnr.Values,'k--.','MarkerSize',10); hold on
plot(rnr.BinEdges(2:end), rnr.Values,'b--.','MarkerSize',10); hold on
plot(nrr.BinEdges(2:end), nrr.Values,'c--.','MarkerSize',10);

legend('R-R','NR-NR','R-NR','NR-R');
title('Population Coupling Distributions', 'FontSize',16)
xlabel('Population Coupling [au]');
ylabel('Probability Density');
set(gca, 'FontSize',16);
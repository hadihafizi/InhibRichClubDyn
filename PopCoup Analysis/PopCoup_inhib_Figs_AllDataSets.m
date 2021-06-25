
clear all

dataSets = {'02-0','03-0','03-1','04-0','04-1','05-0','05-1','06-0',...
    '06-1','07-0','08-0','08-1','09-0','09-1'};
% dataSets = {'139', '151', '152', '163', '165', '168', '174'};
EIdatasets = [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23, 24, 25];
Coup_r_r_norm = [];
Coup_r_e_norm = [];
Coup_r_i_norm = [];
Coup_nr_nr_norm = [];
Coup_nr_e_norm = [];
Coup_nr_i_norm = [];
Coup_r_nr_norm = [];
Coup_rnr_e_norm = [];
Coup_rnr_i_norm = [];
Coup_nr_r_norm = [];
Coup_nrr_e_norm = [];
Coup_nrr_i_norm = [];

for i = 1:length(dataSets)

    % E/I and RC indices
    load([dataSets{i},'\PDF_NoSpur_thr45_',dataSets{i},'.mat'])
    load([dataSets{i},'\wgts_1_16ms.mat'])
    W = PDF.*wgt;oute = sum(W,2);inte = sum(W,1);tote = inte + oute';
    [A, B] = sort(tote,'descend');
    num_rich = ceil(0.2*length(A));
    Rich_idx = B(1:num_rich);
    nonRich_idx = B(num_rich+1:end);

    exc_thresh = .9;
    inh_thresh = .9;
    %in-vitro
    load(['C:\Users\Mahzad\Box Sync\Research\Beggs Lab\Projects\InhibID\Ian Stevenson lab\Neuron Type Analysis\Neuron Type Analysis\in vitro\DataSet', num2str(EIdatasets(i)), '.mat']);
    %in-vivo
%     load(['C:\Users\Mahzad\Box Sync\Research\Beggs Lab\Projects\InhibID\Ian Stevenson lab\Neuron Type Analysis\Neuron Type Analysis\in vivo\unitactivity_Mouse ', num2str(dataSets{i}), '.mat']);
    excprobs = extractfield(neurontype, 'excprob');
    inhprobs = extractfield(neurontype,'inhprob');
    % excConRatio = extractfield(neurontype, 'ExcConRatio');
    ExcIdx = find(excprobs > exc_thresh);
    InhIdx = find(inhprobs > inh_thresh);
    neuron_ids = 1:numel(neurontype);
    unidentified = neuron_ids(~ismember(neuron_ids, ExcIdx) & ~ismember(neuron_ids, InhIdx));
    
    R_E = find(ismember(Rich_idx, intersect(Rich_idx, ExcIdx, 'stable')));
    R_I = find(ismember(Rich_idx, intersect(Rich_idx, InhIdx, 'stable')));
    NR_E = find(ismember(nonRich_idx, intersect(nonRich_idx, ExcIdx, 'stable')));
    NR_I = find(ismember(nonRich_idx, intersect(nonRich_idx, InhIdx, 'stable')));

    % slice PopCoup values to E/I vs RC
    % in-vitro
    load(['PopCoupAll_', dataSets{i}, '_cor_toteRC.mat'])
    %in-vivo
%     load(['PopCoupAll_', dataSets{i}, '_cor_toteRC_bs3ms.mat'])
    [Coup_r_r, Coup_r_nr, Coup_nr_r, Coup_nr_nr] = PopCoup_Figs(Coup_r_r, Coup_r_nr, Coup_nr_r, Coup_nr_nr, dataSets{i}, 0);
    Coup_r_r_norm = [Coup_r_r_norm;Coup_r_r];
    Coup_r_e_norm = [Coup_r_e_norm;Coup_r_r(R_E,:)];
    Coup_r_i_norm = [Coup_r_i_norm;Coup_r_r(R_I,:)];
    Coup_nr_nr_norm = [Coup_nr_nr_norm;Coup_nr_nr];
    Coup_nr_e_norm = [Coup_nr_e_norm;Coup_nr_nr(NR_E,:)];
    Coup_nr_i_norm = [Coup_nr_i_norm;Coup_nr_nr(NR_I,:)];
    Coup_r_nr_norm = [Coup_r_nr_norm;Coup_r_nr];
    Coup_rnr_e_norm = [Coup_rnr_e_norm;Coup_r_nr(R_E,:)];
    Coup_rnr_i_norm = [Coup_rnr_i_norm;Coup_r_nr(R_I,:)];
    Coup_nr_r_norm = [Coup_nr_r_norm;Coup_nr_r];
    Coup_nrr_e_norm = [Coup_nrr_e_norm;Coup_nr_r(NR_E,:)];
    Coup_nrr_i_norm = [Coup_nrr_i_norm;Coup_nr_r(NR_I,:)];
end

%% Actual Data
nbin = 100;
maxCoup = max([max(Coup_r_r_norm(:,:)),max(Coup_r_nr_norm(:,:)),...
    max(Coup_nr_nr_norm(:,:)),max(Coup_nr_r_norm(:,:))]);
minCoup = min([min(Coup_r_r_norm(:,:)),min(Coup_r_nr_norm(:,:)),...
    min(Coup_nr_nr_norm(:,:)),min(Coup_nr_r_norm(:,:))]);
bnedg = minCoup:((maxCoup-minCoup)/nbin):maxCoup;
%%
[rrcount, rrcenter] = hist(Coup_r_r_norm(:,51),bnedg); 
[recount, recenter] = hist(Coup_r_e_norm(:,51),bnedg); 
[ricount, ricenter] = hist(Coup_r_i_norm(:,51),bnedg); 
[nrnrcount, nrnrcenter] = hist(Coup_nr_nr_norm(:,51),bnedg); 
[nrecount, nrecenter] = hist(Coup_nr_e_norm(:,51),bnedg); 
[nricount, nricenter] = hist(Coup_nr_i_norm(:,51),bnedg); 

rrcount = rrcount/sum(rrcount);
recount = recount/sum(recount);
ricount = ricount/sum(ricount);
nrnrcount = nrnrcount/sum(nrnrcount);
nrecount = nrecount/sum(nrecount);
nricount = nricount/sum(nricount);

figure;
plot(rrcenter,rrcount,'r--.','MarkerSize',10); hold on
plot(recenter,recount,'g--.','MarkerSize',10); hold on
plot(ricenter,ricount,'m--.','MarkerSize',10); hold on
plot(nrnrcenter,nrnrcount,'b--.','MarkerSize',10); hold on
plot(nrecenter,nrecount,'c--.','MarkerSize',10); hold on
plot(nricenter,nricount,'k--.','MarkerSize',10); hold on

legend('R-R','RE-R','RI-R','NR-NR','NRE-NR','NRI-NR');
title('Population Coupling Distributions', 'FontSize',16)
xlabel('Population Coupling [au]');
ylabel('Probability');
set(gca, 'FontSize',16);

%% R-R and NR-NR
ax = figure(); hold on
rr = histogram(Coup_r_r_norm(:,51),bnedg, 'normalization', 'pdf', 'DisplayStyle','stairs', 'LineWidth', 2);
re = histogram(Coup_r_e_norm(:,51),bnedg, 'normalization', 'pdf', 'DisplayStyle','stairs', 'LineWidth', 2);
ri = histogram(Coup_r_i_norm(:,51),bnedg, 'normalization', 'pdf', 'DisplayStyle','stairs', 'LineWidth', 2);
nrnr = histogram(Coup_nr_nr_norm(:,51),bnedg, 'normalization', 'pdf', 'DisplayStyle','stairs', 'LineWidth', 2);
nre = histogram(Coup_nr_e_norm(:,51),bnedg, 'normalization', 'pdf', 'DisplayStyle','stairs', 'LineWidth', 2);
nri = histogram(Coup_nr_i_norm(:,51),bnedg, 'normalization', 'pdf', 'DisplayStyle','stairs', 'LineWidth', 2);

legend('R-R','RE-R','RI-R','NR-NR','NRE-NR','NRI-NR');
title('Population Coupling Distributions', 'FontSize',16)
xlabel('Population Coupling [au]');
ylabel('Probability Density');
set(gca, 'FontSize',16);

figure;
plot(rr.BinEdges(2:end), rr.Values,'r--.','MarkerSize',10); hold on
plot(re.BinEdges(2:end), re.Values,'g--.','MarkerSize',10); hold on
plot(ri.BinEdges(2:end), ri.Values,'m--.','MarkerSize',10);
plot(nrnr.BinEdges(2:end), nrnr.Values,'b--.','MarkerSize',10); hold on
plot(nre.BinEdges(2:end), nre.Values,'c--.','MarkerSize',10); hold on
plot(nri.BinEdges(2:end), nri.Values,'k--.','MarkerSize',10); hold on

legend('R-R','RE-R','RI-R','NR-NR','NRE-NR','NRI-NR');
title('Population Coupling Distributions', 'FontSize',16)
xlabel('Population Coupling [au]');
ylabel('Probability Density');
set(gca, 'FontSize',16);

%% R-NR and NR-R
ax = figure(); hold on
rnr = histogram(Coup_r_nr_norm(:,51),bnedg, 'normalization', 'pdf', 'DisplayStyle','stairs', 'LineWidth', 2);
rnre = histogram(Coup_rnr_e_norm(:,51),bnedg, 'normalization', 'pdf', 'DisplayStyle','stairs', 'LineWidth', 2);
rnri = histogram(Coup_rnr_i_norm(:,51),bnedg, 'normalization', 'pdf', 'DisplayStyle','stairs', 'LineWidth', 2);
nrr = histogram(Coup_nr_r_norm(:,51),bnedg, 'normalization', 'pdf', 'DisplayStyle','stairs', 'LineWidth', 2);
nrre = histogram(Coup_nrr_e_norm(:,51),bnedg, 'normalization', 'pdf', 'DisplayStyle','stairs', 'LineWidth', 2);
nrri = histogram(Coup_nrr_i_norm(:,51),bnedg, 'normalization', 'pdf', 'DisplayStyle','stairs', 'LineWidth', 2);

legend('R-NR','RE-NR','RI-NR','NR-R','NRE-R','NRI-R');
title('Population Coupling Distributions', 'FontSize',16)
xlabel('Population Coupling [au]');
ylabel('Probability Density');
set(gca, 'FontSize',16);

figure;
plot(rnr.BinEdges(2:end), rnr.Values,'r--.','MarkerSize',10); hold on
plot(rnre.BinEdges(2:end), rnre.Values,'g--.','MarkerSize',10); hold on
plot(rnri.BinEdges(2:end), rnri.Values,'m--.','MarkerSize',10);
plot(nrr.BinEdges(2:end), nrr.Values,'b--.','MarkerSize',10); hold on
plot(nrre.BinEdges(2:end), nrre.Values,'c--.','MarkerSize',10); hold on
plot(nrri.BinEdges(2:end), nrri.Values,'k--.','MarkerSize',10); hold on

legend('R-NR','RE-NR','RI-NR','NR-R','NRE-R','NRI-R');
title('Population Coupling Distributions', 'FontSize',16)
xlabel('Population Coupling [au]');
ylabel('Probability Density');
set(gca, 'FontSize',16);

%% 

[A1, B1] = size(Coup_r_r_norm);
interval = (B1-1)/2;
TLag = -interval:interval;

AUC_Pre_r_r = sum((Coup_r_r_norm(:,1:(length(TLag)-1)/2)),2);
AUC_Post_r_r = sum((Coup_r_r_norm(:,1+(length(TLag)-1)/2:end)),2);
driver_r = (AUC_Pre_r_r - AUC_Post_r_r)./(AUC_Pre_r_r + AUC_Post_r_r);
AUC_Pre_r_e = sum((Coup_r_e_norm(:,1:(length(TLag)-1)/2)),2);
AUC_Post_r_e = sum((Coup_r_e_norm(:,1+(length(TLag)-1)/2:end)),2);
driver_r_e = (AUC_Pre_r_e - AUC_Post_r_e)./(AUC_Pre_r_e + AUC_Post_r_e);
AUC_Pre_r_i = sum((Coup_r_i_norm(:,1:(length(TLag)-1)/2)),2);
AUC_Post_r_i = sum((Coup_r_i_norm(:,1+(length(TLag)-1)/2:end)),2);
driver_r_i = (AUC_Pre_r_i - AUC_Post_r_i)./(AUC_Pre_r_i + AUC_Post_r_i);
AUC_Pre_nr_nr = sum((Coup_nr_nr_norm(:,1:(length(TLag)-1)/2)),2);
AUC_Post_nr_nr = sum((Coup_nr_nr_norm(:,1+(length(TLag)-1)/2:end)),2);
driver_nr = (AUC_Pre_nr_nr - AUC_Post_nr_nr)./(AUC_Pre_nr_nr + AUC_Post_nr_nr);
AUC_Pre_nr_e = sum((Coup_nr_e_norm(:,1:(length(TLag)-1)/2)),2);
AUC_Post_nr_e = sum((Coup_nr_e_norm(:,1+(length(TLag)-1)/2:end)),2);
driver_nr_e = (AUC_Pre_nr_e - AUC_Post_nr_e)./(AUC_Pre_nr_e + AUC_Post_nr_e);
AUC_Pre_nr_i = sum((Coup_nr_i_norm(:,1:(length(TLag)-1)/2)),2);
AUC_Post_nr_i = sum((Coup_nr_i_norm(:,1+(length(TLag)-1)/2:end)),2);
driver_nr_i = (AUC_Pre_nr_i - AUC_Post_nr_i)./(AUC_Pre_nr_i + AUC_Post_nr_i);
AUC_Pre_r_nr = sum((Coup_r_nr_norm(:,1:(length(TLag)-1)/2)),2);
AUC_Post_r_nr = sum((Coup_r_nr_norm(:,1+(length(TLag)-1)/2:end)),2);
driver_rnr = (AUC_Pre_r_nr - AUC_Post_r_nr)./(AUC_Pre_r_nr + AUC_Post_r_nr);
AUC_Pre_rnr_e = sum((Coup_rnr_e_norm(:,1:(length(TLag)-1)/2)),2);
AUC_Post_rnr_e = sum((Coup_rnr_e_norm(:,1+(length(TLag)-1)/2:end)),2);
driver_rnr_e = (AUC_Pre_rnr_e - AUC_Post_rnr_e)./(AUC_Pre_rnr_e + AUC_Post_rnr_e);
AUC_Pre_rnr_i = sum((Coup_rnr_i_norm(:,1:(length(TLag)-1)/2)),2);
AUC_Post_rnr_i = sum((Coup_rnr_i_norm(:,1+(length(TLag)-1)/2:end)),2);
driver_rnr_i = (AUC_Pre_rnr_i - AUC_Post_rnr_i)./(AUC_Pre_rnr_i + AUC_Post_rnr_i);
AUC_Pre_nr_r = sum((Coup_nr_r_norm(:,1:(length(TLag)-1)/2)),2);
AUC_Post_nr_r = sum((Coup_nr_r_norm(:,1+(length(TLag)-1)/2:end)),2);
driver_nrr = (AUC_Pre_nr_r - AUC_Post_nr_r)./(AUC_Pre_nr_r + AUC_Post_nr_r);
AUC_Pre_nrr_e = sum((Coup_nrr_e_norm(:,1:(length(TLag)-1)/2)),2);
AUC_Post_nrr_e = sum((Coup_nrr_e_norm(:,1+(length(TLag)-1)/2:end)),2);
driver_nrr_e = (AUC_Pre_nrr_e - AUC_Post_nrr_e)./(AUC_Pre_nrr_e + AUC_Post_nrr_e);
AUC_Pre_nrr_i = sum((Coup_nrr_i_norm(:,1:(length(TLag)-1)/2)),2);
AUC_Post_nrr_i = sum((Coup_nrr_i_norm(:,1+(length(TLag)-1)/2:end)),2);
driver_nrr_i = (AUC_Pre_nrr_i - AUC_Post_nrr_i)./(AUC_Pre_nrr_i + AUC_Post_nrr_i);

figure;
scatter(AUC_Pre_nr_e, AUC_Post_nr_e, 'c'); hold on
scatter(AUC_Pre_nr_i, AUC_Post_nr_i, 'k'); 
% scatter(AUC_Pre_nrr_e, AUC_Post_nrr_e, 'r');
% scatter(AUC_Pre_nrr_i, AUC_Post_nrr_i, 'rx');
% scatter(AUC_Pre_rnr_i, AUC_Post_rnr_i, 'y');
% scatter(AUC_Pre_rnr_e, AUC_Post_rnr_e, 'b');
scatter(AUC_Pre_r_i, AUC_Post_r_i, 'm');
scatter(AUC_Pre_r_e, AUC_Post_r_e, 'g');
xlabel('AUC_{pre}')
ylabel('AUC_{post}')
legend('NR_e','NR_i','R_i','R_e')
grid on
figure;
scatter(driver_nr_e, AUC_Post_nr_e, 'c'); hold on
scatter(driver_nr_i, AUC_Post_nr_i, 'k'); 
% scatter(driver_nrr_e, AUC_Post_nrr_e, 'r');
% scatter(driver_nrr_i, AUC_Post_nrr_i, 'rx');
% scatter(driver_rnr_i, AUC_Post_rnr_i, 'y');
% scatter(driver_rnr_e, AUC_Post_rnr_e, 'b');
scatter(driver_r_i, AUC_Post_r_i, 'm');
scatter(driver_r_e, AUC_Post_r_e, 'g');
xlabel('(AUC_{pre}-AUC_{post})/(AUC_{pre}+AUC_{post})')
ylabel('AUC_{post}')
legend('NR_e','NR_i','R_i','R_e')
grid on

figure; hold on
histogram(driver_r_e, 'Normalization', 'pdf');
histogram(driver_r_i, 'Normalization', 'pdf');
histogram(driver_nr_e, 'Normalization', 'pdf');
histogram(driver_nr_i, 'Normalization', 'pdf');
% histogram(driver_rnr_e, 'Normalization', 'pdf');
% histogram(driver_rnr_i, 'Normalization', 'pdf');
% histogram(driver_nrr_e, 'Normalization', 'pdf');
% histogram(driver_nrr_i, 'Normalization', 'pdf');
xlabel('(AUC_{pre}-AUC_{post})/(AUC_{pre}+AUC_{post})')
ylabel('pdf')
legend('R_e','R_i','NR_e','NR_i')







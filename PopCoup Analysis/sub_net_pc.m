% qwe = Coup(1,:);
% ii = 1;kk=1;
% while ii+20 < length(qwe)
%     spk(kk) = sum(qwe(ii:ii+20));
%     ii = ii+20;
%     kk = kk+1;
% end
% plot(spk)
clear all
load Coup_r_r.mat;load Coup_nr_nr.mat;
load Coup_r_nr.mat;load Coup_nr_r.mat;
[A1 B1] = size(Coup_r_r);[A2 B2] = size(Coup_nr_nr);

Coup_r_r_shuf_sum = zeros(A1,B1);
Coup_nr_nr_shuf_sum = zeros(A2,B2);
Coup_r_nr_shuf_sum = zeros(A1,B1);
Coup_nr_r_shuf_sum = zeros(A2,B2);

for ii=1:20
    str = ['cd shuf',num2str(ii)]; eval(str);
    str = ['load Coup_r_r_shuf',num2str(ii),'.mat']; eval(str);
    str = ['load Coup_nr_nr_shuf',num2str(ii),'.mat']; eval(str);
    str = ['load Coup_r_nr_shuf',num2str(ii),'.mat']; eval(str);
    str = ['load Coup_nr_r_shuf',num2str(ii),'.mat']; eval(str);
    
    Coup_r_r_shuf_sum = Coup_r_r_shuf + Coup_r_r_shuf_sum;
    Coup_nr_nr_shuf_sum = Coup_nr_nr_shuf + Coup_nr_nr_shuf_sum;
    Coup_r_nr_shuf_sum = Coup_r_nr_shuf + Coup_r_nr_shuf_sum;
    Coup_nr_r_shuf_sum = Coup_nr_r_shuf + Coup_nr_r_shuf_sum;
    
    cd ..
end
Coup_r_r_shuf_sum = Coup_r_r_shuf_sum/ii;
Coup_nr_nr_shuf_sum = Coup_nr_nr_shuf_sum/ii;
Coup_r_nr_shuf_sum = Coup_r_nr_shuf_sum/ii;
Coup_nr_r_shuf_sum = Coup_nr_r_shuf_sum/ii;
% cd shuf1
% load Coup_r_r_shuf.mat; load Coup_nr_nr_shuf.mat;
% load Coup_r_nr_shuf.mat; load Coup_nr_r_shuf.mat
% cd ..

% normalize for the number of nodes and subtract baseline. 
for ii = 1:A1
    Coup_r_r_norm(ii,:) = (Coup_r_r(ii,:)-Coup_r_r(ii,1))./(A1-1);
    Coup_r_r_shuf_norm(ii,:) = (Coup_r_r_shuf_sum(ii,:)-Coup_r_r_shuf_sum(ii,1))./(A1-1);
end
TLag = [-500:500];
%figure
% subplot(1,2,1);plot([-500:500],aa')
% subplot(1,2,2);plot([-500:500],aa_shuf')
% [s1,s2] = size(Coup_r_r);
figure(202);
subplot(2,2,1);
surf(TLag,[1:A1],Coup_r_r_norm); colormap('jet');colorbar; shading('interp')
% subplot(1,2,2);surf(aa_shuf,xx); colormap('jet');colorbar; shading('interp')
title('Rich-Rich', 'FontSize',16)
set(gca, 'FontSize',16);
xlabel('Time-lag [ms]')
ylabel('Neuron ID')
view(2)
% savefig('pop_coup_r_r.fig');
figure;
histogram(Coup_r_r_norm(:,501),20,'normalization','probability'); hold on
histogram(Coup_r_r_shuf_norm(:,501),'normalization','probability'); hold on
legend('RR-Actual','RR-Shuff');
title('R-R Act-shuf', 'FontSize',16)
set(gca, 'FontSize',16);
% savefig('pop_coup_hist_rrActVSrrShuff.fig');


for ii = 1:A2
    Coup_nr_nr_norm(ii,:) = (Coup_nr_nr(ii,:)-Coup_nr_nr(ii,1))./(A2-1);
    Coup_nr_nr_shuf_norm(ii,:) = (Coup_nr_nr_shuf_sum(ii,:)-Coup_nr_nr_shuf_sum(ii,1))./(A2-1);
end
% figure
% subplot(1,2,1);plot([-500:500],bb')
% subplot(1,2,2);plot([-500:500],bb_shuf')
% [s1,s2] = size(Coup_nr_nr);
figure(202);
subplot(2,2,2);
surf(TLag,[1:A2],Coup_nr_nr_norm); colormap('jet');h = colorbar; shading('interp')
% lims = get(h,'Limits'); clear h
%subplot(1,2,2);surf(bb_shuf); colormap('jet');h1 = colorbar; shading('interp')
% set(h1,'Limits',[lims])
title('nonRich-nonRich', 'FontSize',16)
set(gca, 'FontSize',16);
xlabel('Time-lag [ms]')
ylabel('Neuron ID')
view(2)
% savefig('pop_coup_nr_nr.fig');
figure;
histogram(Coup_r_r_norm(:,501),20,'normalization','probability'); hold on
histogram(Coup_nr_nr_norm(:,501),'normalization','probability'); hold on
legend('R-R','NR-NR');
title('R-R vs NR-NR', 'FontSize',16)
set(gca, 'FontSize',16);
% savefig('pop_coup_hist_rrVSnrnr.fig');

for ii = 1:A1
    Coup_r_nr_norm(ii,:) = (Coup_r_nr(ii,:)-Coup_r_nr(ii,1))./(A2-1);
    Coup_r_nr_shuf_norm(ii,:) = (Coup_r_nr_shuf_sum(ii,:)-Coup_r_nr_shuf_sum(ii,1))./(A2-1);
end
% figure
% subplot(1,2,1);plot([-500:500],cc')
% subplot(1,2,2);plot([-500:500],cc_shuf')
figure(202);
% [s1,s2] = size(Coup_r_nr);
subplot(2,2,3);
surf(TLag,[1:A1],Coup_r_nr_norm); colormap('jet');colorbar; shading('interp')
%subplot(1,2,2);surf(cc_shuf); colormap('jet');colorbar; shading('interp')
title('Rich-nonRich', 'FontSize',16)
set(gca, 'FontSize',16);
xlabel('Time-lag [ms]')
ylabel('Neuron ID')
view(2)
% savefig('pop_coup_r_nr.fig');
figure;
histogram(Coup_r_r_norm(:,501),20,'normalization','probability'); hold on
histogram(Coup_r_nr_norm(:,501),'normalization','probability'); hold on
legend('R-R','R-NR');
title('R-R vs R-NR', 'FontSize',16)
set(gca, 'FontSize',16);
% savefig('pop_coup_hist_rrVSrnr.fig');

for ii = 1:A2
    Coup_nr_r_norm(ii,:) = (Coup_nr_r(ii,:)-Coup_nr_r(ii,1))./(A1-1);
    Coup_nr_r_shuf_norm(ii,:) = (Coup_nr_r_shuf_sum(ii,:)-Coup_nr_r_shuf_sum(ii,1))./(A1-1);
end
% figure
% subplot(1,2,1);plot([-500:500],Coup_nr_r_norm')
% subplot(1,2,2);plot([-500:500],Coup_nr_r_shuf_norm')
figure(202);
subplot(2,2,4);
surf(TLag,[1:A2],Coup_nr_r_norm); colormap('jet');colorbar; shading('interp'); %view(2)
set(gca, 'FontSize',16);
% subplot(1,2,2);
% surf(TLag,[1:A2],Coup_nr_r_shuf_norm); colormap('jet');colorbar; shading('interp')
title('nonRich-Rich', 'FontSize',16)
xlabel('Time-lag [ms]')
ylabel('Neuron ID')
view(2)
% savefig('pop_coup_nr_r.fig');
figure;
histogram(Coup_nr_nr_norm(:,501),'normalization','probability'); hold on
histogram(Coup_nr_r_norm(:,501),20,'normalization','probability'); hold on
legend('NR-NR','NR-R')
title('NR-NR vs NR-R', 'FontSize',16)
set(gca, 'FontSize',16);
% savefig('pop_coup_hist_nrnrVSnrr.fig');


F1 = figure;
NumBin = 40;
h1 = histogram(Coup_nr_r_norm(:,501),NumBin,'normalization','probability'); hold on
bnedg = h1.BinEdges;hold on
histogram(Coup_r_r_norm(:,501),bnedg,'normalization','probability'); hold on
legend('NR-R','R-R')
title('NR-R vs R-R', 'FontSize',16)
set(gca, 'FontSize',16);
% cd Figures
% savefig('pop_coup_hist_nrrVSrr.fig');
% saveas(F1,['pop_coup_hist_nrrVSrr.jpg'])
% cd ..

%% PopCoupAll

nbin = 50;
maxCoup = max([max(Coup_r_r_norm(:,:)),max(Coup_nr_r_norm(:,:)),...
    max(Coup_nr_nr_norm(:,:)),max(Coup_r_nr_norm(:,:))]);
minCoup = min([min(Coup_r_r_norm(:,:)),min(Coup_nr_r_norm(:,:)),...
    min(Coup_nr_nr_norm(:,:)),min(Coup_r_nr_norm(:,:))]);
bnedg = minCoup:((maxCoup-minCoup)/nbin):maxCoup;

[rrcount, rrcenter] = hist(Coup_r_r_norm(:,501),bnedg); 
[nrrcount, nrrcenter] = hist(Coup_nr_r_norm(:,501),bnedg); 
[nrnrcount, nrnrcenter] = hist(Coup_nr_nr_norm(:,501),bnedg); 
[rnrcount, rnrcenter] = hist(Coup_r_nr_norm(:,501),bnedg); 

rrcount = rrcount/sum(rrcount);
nrnrcount = nrnrcount/sum(nrnrcount);
rnrcount = rnrcount/sum(rnrcount);
nrrcount = nrrcount/sum(nrrcount);

frr = fit(rrcenter.',rrcount.','gauss1');
fnrnr = fit(nrnrcenter.',nrnrcount.','gauss1');
frnr = fit(rnrcenter.',rnrcount.','gauss1');
fnrr = fit(nrrcenter.',nrrcount.','gauss2');

figure;
plot(frr,rrcenter,rrcount); hold on
plot(fnrnr,nrnrcenter,nrnrcount); hold on
plot(frnr,rnrcenter,rnrcount); hold on
plot(fnrr,nrrcenter,nrrcount); hold on

legend('','R-R','','NR-NR','','R-NR','','NR-R');
title('Population Coupling Distributions', 'FontSize',16)
xlabel('Population Coupling [au]');
ylabel('Probability');
set(gca, 'FontSize',16);



%% Driver/Follower (Asymmetry of pop_coup curves)

AUC_Pre_nr_r = sum((Coup_nr_r_norm(:,1:(length(TLag)-1)/2)),2);
AUC_Post_nr_r = sum((Coup_nr_r_norm(:,1+(length(TLag)-1)/2:end)),2);

% for ii=1:20
    AUC_Pre_nr_r_shuf = sum((Coup_nr_r_shuf_norm(:,1:(length(TLag)-1)/2)),2);
    AUC_Post_nr_r_shuf = sum((Coup_nr_r_shuf_norm(:,1+(length(TLag)-1)/2:end)),2);
% end

figure;scatter(AUC_Pre_nr_r,AUC_Post_nr_r,'o', 'LineWidth',2);hold on
scatter(AUC_Pre_nr_r_shuf,AUC_Post_nr_r_shuf,'o', 'LineWidth',2);hold on
plot(AUC_Pre_nr_r,AUC_Pre_nr_r, 'LineWidth',2)
legend('Actual','Shuffled')
axis([min(min(AUC_Pre_nr_r),min(AUC_Post_nr_r)) max(max(AUC_Pre_nr_r),max(AUC_Post_nr_r)) min(min(AUC_Pre_nr_r),min(AUC_Post_nr_r)) max(max(AUC_Pre_nr_r),max(AUC_Post_nr_r))])
title('Asymmetry in Population Couplings (NR-R)','FontSize',16);xlabel('Pre-Spike','FontSize',16);ylabel('Post-Spike','FontSize',16);
set(gca, 'FontSize',16);
% savefig('pop_coup_asym_nr_r.fig');
%%%
AUC_Pre_nr_nr = sum((Coup_nr_nr_norm(:,1:(length(TLag)-1)/2)),2);
AUC_Post_nr_nr = sum((Coup_nr_nr_norm(:,1+(length(TLag)-1)/2:end)),2);

AUC_Pre_nr_nr_shuf = sum((Coup_nr_nr_shuf_norm(:,1:(length(TLag)-1)/2)),2);
AUC_Post_nr_nr_shuf = sum((Coup_nr_nr_shuf_norm(:,1+(length(TLag)-1)/2:end)),2);

figure;scatter(AUC_Pre_nr_nr,AUC_Post_nr_nr,'o', 'LineWidth',2);hold on
scatter(AUC_Pre_nr_nr_shuf,AUC_Post_nr_nr_shuf,'o', 'LineWidth',2);hold on
plot(AUC_Pre_nr_nr,AUC_Pre_nr_nr, 'LineWidth',2)
legend('Actual','Shuffled')
axis([min(min(AUC_Pre_nr_nr),min(AUC_Post_nr_nr)) max(max(AUC_Pre_nr_nr),max(AUC_Post_nr_nr)) min(min(AUC_Pre_nr_nr),min(AUC_Post_nr_nr)) max(max(AUC_Pre_nr_nr),max(AUC_Post_nr_nr))])
title('Asymmetry in Population Couplings (NR-NR)','FontSize',16);xlabel('Pre-Spike','FontSize',16);ylabel('Post-Spike','FontSize',16);
set(gca, 'FontSize',16);
% savefig('pop_coup_asym_nr_nr.fig');
%%%
AUC_Pre_r_r = sum((Coup_r_r_norm(:,1:(length(TLag)-1)/2)),2);
AUC_Post_r_r = sum((Coup_r_r_norm(:,1+(length(TLag)-1)/2:end)),2);

AUC_Pre_r_r_shuf = sum((Coup_r_r_shuf_norm(:,1:(length(TLag)-1)/2)),2);
AUC_Post_r_r_shuf = sum((Coup_r_r_shuf_norm(:,1+(length(TLag)-1)/2:end)),2);

figure;scatter(AUC_Pre_r_r,AUC_Post_r_r,'o', 'LineWidth',2);hold on
scatter(AUC_Pre_r_r_shuf,AUC_Post_r_r_shuf,'o', 'LineWidth',2);hold on
plot(AUC_Pre_r_r,AUC_Pre_r_r, 'LineWidth',2)
legend('Actual','Shuffled')
axis([min(min(AUC_Pre_r_r),min(AUC_Post_r_r)) max(max(AUC_Pre_r_r),max(AUC_Post_r_r)) min(min(AUC_Pre_r_r),min(AUC_Post_r_r)) max(max(AUC_Pre_r_r),max(AUC_Post_r_r))])
title('Asymmetry in Population Couplings (R-R)','FontSize',16);xlabel('Pre-Spike','FontSize',16);ylabel('Post-Spike','FontSize',16);
set(gca, 'FontSize',16);
% savefig('pop_coup_asym_r_r.fig');
%%%
AUC_Pre_r_nr = sum((Coup_r_nr_norm(:,1:(length(TLag)-1)/2)),2);
AUC_Post_r_nr = sum((Coup_r_nr_norm(:,1+(length(TLag)-1)/2:end)),2);

AUC_Pre_r_nr_shuf = sum((Coup_r_nr_shuf_norm(:,1:(length(TLag)-1)/2)),2);
AUC_Post_r_nr_shuf = sum((Coup_r_nr_shuf_norm(:,1+(length(TLag)-1)/2:end)),2);

figure;scatter(AUC_Pre_r_nr,AUC_Post_r_nr,'o', 'LineWidth',2);hold on
scatter(AUC_Pre_r_nr_shuf,AUC_Post_r_nr_shuf,'o', 'LineWidth',2);hold on
plot(AUC_Pre_r_nr,AUC_Pre_r_nr, 'LineWidth',2)
legend('Actual','Shuffled')
axis([min(min(AUC_Pre_r_nr),min(AUC_Post_r_nr)) max(max(AUC_Pre_r_nr),max(AUC_Post_r_nr)) min(min(AUC_Pre_r_nr),min(AUC_Post_r_nr)) max(max(AUC_Pre_r_nr),max(AUC_Post_r_nr))])
title('Asymmetry in Population Couplings (R-NR)','FontSize',16);xlabel('Pre-Spike','FontSize',16);ylabel('Post-Spike','FontSize',16);
set(gca, 'FontSize',16);
% savefig('pop_coup_asym_r_nr.fig');





% function [Coup_r_r_norm, Coup_r_nr_norm, Coup_nr_r_norm, Coup_nr_nr_norm, ...
%     Coup_r_r_shuf_norm, Coup_r_nr_shuf_norm, Coup_nr_r_shuf_norm, Coup_nr_nr_shuf_norm] = ...
%     PopCoup_Figs(Coup_r_r,Coup_r_nr,Coup_nr_r,Coup_nr_nr, dataSet,varargin)
function [Coup_r_r_norm, Coup_r_nr_norm, Coup_nr_r_norm, Coup_nr_nr_norm] = ...
    PopCoup_Figs(Coup_r_r,Coup_r_nr,Coup_nr_r,Coup_nr_nr, dataSet,varargin)
% varargin:
%         - Figs [0/1]: whether or not to plot the figures.
% varargout:
%         - Coup_r_r_norm  
%         - Coup_r_nr_norm  
%         - Coup_nr_r_norm  
%         - Coup_nr_nr_norm  
% clear all
% dataSet = '05-1';
% suffix = 'rwrandpermSame1st';  %'rw' for rewired nets; '' for originals; '_cor' for
% load([dataSet,'/PopCoupAll_',dataSet,suffix,'Simk',PFeigStr,'.mat'])
[A1, B1] = size(Coup_r_r);[A2, B2] = size(Coup_nr_nr);
shuf = 0;
interval = (B1-1)/2;
TLag = -interval:interval;

% normalize for the number of nodes and subtract baseline.
for ii = 1:A1
%     Coup_r_r_norm(ii,:) = (Coup_r_r(ii,:)-Coup_r_r(ii,1))./(A1-1);
%     Coup_r_nr_norm(ii,:) = (Coup_r_nr(ii,:)-Coup_r_nr(ii,1))./(A2-1);
    Coup_r_r_norm(ii,:) = (Coup_r_r(ii,:))./(A1-1);
    Coup_r_nr_norm(ii,:) = (Coup_r_nr(ii,:))./(A2-1);
end
for ii = 1:A2
%     Coup_nr_nr_norm(ii,:) = (Coup_nr_nr(ii,:)-Coup_nr_nr(ii,1))./(A2-1);
%     Coup_nr_r_norm(ii,:) = (Coup_nr_r(ii,:)-Coup_nr_r(ii,1))./(A1-1);
    Coup_nr_nr_norm(ii,:) = (Coup_nr_nr(ii,:))./(A2-1);
    Coup_nr_r_norm(ii,:) = (Coup_nr_r(ii,:))./(A1-1);
end

if shuf
    Coup_r_r_shuf_sum = zeros(A1,B1);
    Coup_nr_nr_shuf_sum = zeros(A2,B2);
    Coup_r_nr_shuf_sum = zeros(A1,B1);
    Coup_nr_r_shuf_sum = zeros(A2,B2);
    
    for ii=1:20
        load(['PopCoupAll_', dataSet, '_cor_toteRC_bs1ms_shuf', num2str(ii), '.mat'])
%         str = ['cd shuf',num2str(ii)]; eval(str);
%         str = ['load Coup_r_r_shuf',num2str(ii),'.mat']; eval(str);
%         str = ['load Coup_nr_nr_shuf',num2str(ii),'.mat']; eval(str);
%         str = ['load Coup_r_nr_shuf',num2str(ii),'.mat']; eval(str);
%         str = ['load Coup_nr_r_shuf',num2str(ii),'.mat']; eval(str);
        
        Coup_r_r_shuf_sum = Coup_r_r + Coup_r_r_shuf_sum;
        Coup_nr_nr_shuf_sum = Coup_nr_nr + Coup_nr_nr_shuf_sum;
        Coup_r_nr_shuf_sum = Coup_r_nr + Coup_r_nr_shuf_sum;
        Coup_nr_r_shuf_sum = Coup_nr_r + Coup_nr_r_shuf_sum;
        
%         cd ..
    end
    Coup_r_r_shuf_sum = Coup_r_r_shuf_sum/ii;
    Coup_nr_nr_shuf_sum = Coup_nr_nr_shuf_sum/ii;
    Coup_r_nr_shuf_sum = Coup_r_nr_shuf_sum/ii;
    Coup_nr_r_shuf_sum = Coup_nr_r_shuf_sum/ii;

    % normalize for the number of nodes and subtract baseline.
    for ii = 1:A1
        Coup_r_r_shuf_norm(ii,:) = (Coup_r_r_shuf_sum(ii,:)-Coup_r_r_shuf_sum(ii,1))./(A1-1);
        Coup_r_nr_shuf_norm(ii,:) = (Coup_r_nr_shuf_sum(ii,:)-Coup_r_nr_shuf_sum(ii,1))./(A2-1);
    end
    for ii = 1:A2
        Coup_nr_nr_shuf_norm(ii,:) = (Coup_nr_nr_shuf_sum(ii,:)-Coup_nr_nr_shuf_sum(ii,1))./(A2-1);
        Coup_nr_r_shuf_norm(ii,:) = (Coup_nr_r_shuf_sum(ii,:)-Coup_nr_r_shuf_sum(ii,1))./(A1-1);
    end
end

% if nargout > 0
%     varargout = {Coup_r_r_norm, Coup_r_nr_norm, ...
%         Coup_nr_r_norm, Coup_nr_nr_norm};
% end

%% Plot Figures
if max(nargin,1) > 5
    Figs = varargin{1};
else
    Figs = 0;
end
if Figs
    %{
figure(202);
subplot(2,2,1);
surf(TLag,[1:A1],Coup_r_r_norm); colormap('jet');colorbar; shading('interp')
title('Rich-Rich', 'FontSize',16)
set(gca, 'FontSize',16);
xlabel('Time-lag [ms]')
ylabel('Neuron ID')
view(2)

subplot(2,2,2);
surf(TLag,[1:A2],Coup_nr_nr_norm); colormap('jet');h = colorbar; shading('interp')
title('nonRich-nonRich', 'FontSize',16)
set(gca, 'FontSize',16);
xlabel('Time-lag [ms]')
ylabel('Neuron ID')
view(2)

subplot(2,2,3);
surf(TLag,[1:A1],Coup_r_nr_norm); colormap('jet');colorbar; shading('interp')
title('Rich-nonRich', 'FontSize',16)
set(gca, 'FontSize',16);
xlabel('Time-lag [ms]')
ylabel('Neuron ID')
view(2)

subplot(2,2,4);
surf(TLag,[1:A2],Coup_nr_r_norm); colormap('jet');colorbar; shading('interp'); %view(2)
set(gca, 'FontSize',16);
title('nonRich-Rich', 'FontSize',16)
xlabel('Time-lag [ms]')
ylabel('Neuron ID')
view(2)
%}
%% PopCoupAll

nbin = 50;
maxCoup = max([max(Coup_r_r_norm(:,:)),max(Coup_nr_r_norm(:,:)),...
    max(Coup_nr_nr_norm(:,:)),max(Coup_r_nr_norm(:,:))]);
minCoup = min([min(Coup_r_r_norm(:,:)),min(Coup_nr_r_norm(:,:)),...
    min(Coup_nr_nr_norm(:,:)),min(Coup_r_nr_norm(:,:))]);
bnedg = minCoup:((maxCoup-minCoup)/nbin):maxCoup;

[rrcount, rrcenter] = hist(Coup_r_r_norm(:,interval+1),bnedg); 
[nrrcount, nrrcenter] = hist(Coup_nr_r_norm(:,interval+1),bnedg); 
[nrnrcount, nrnrcenter] = hist(Coup_nr_nr_norm(:,interval+1),bnedg); 
[rnrcount, rnrcenter] = hist(Coup_r_nr_norm(:,interval+1),bnedg); 

rrcount = rrcount/sum(rrcount);
nrnrcount = nrnrcount/sum(nrnrcount);
rnrcount = rnrcount/sum(rnrcount);
nrrcount = nrrcount/sum(nrrcount);

figure;
%{
frr = fit(rrcenter.',rrcount.','gauss1');
fnrnr = fit(nrnrcenter.',nrnrcount.','gauss1');
frnr = fit(rnrcenter.',rnrcount.','gauss1');
fnrr = fit(nrrcenter.',nrrcount.','gauss2');

figure;
plot(frr,rrcenter,rrcount); hold on
plot(fnrnr,nrnrcenter,nrnrcount); hold on
plot(frnr,rnrcenter,rnrcount); hold on
plot(fnrr,nrrcenter,nrrcount); hold on
%}

plot(rrcenter,rrcount,'r--.','MarkerSize',10); hold on
plot(nrnrcenter,nrnrcount,'k--.','MarkerSize',10); hold on
plot(rnrcenter,rnrcount,'b--.','MarkerSize',10); hold on
plot(nrrcenter,nrrcount,'c--.','MarkerSize',10); hold on

% legend('','R-R','','NR-NR','','R-NR');
legend('R-R','NR-NR','R-NR','NR-R');
title('Population Coupling Distributions', 'FontSize',16)
xlabel('Population Coupling [au]');
ylabel('Probability');
set(gca, 'FontSize',16);
end
%% Driver/Follower (Asymmetry of pop_coup curves)
%{
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

%}



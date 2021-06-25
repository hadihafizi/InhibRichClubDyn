% Calculate the statistical significance of the difference between
% pair-wise similarity of spike trains and their jittered versions. The
% p-values from ks-tests are plotted.
%
% Required files: Sim_bin[X]ms.mat
%                 Sim_jit[X]bin_[X]ms.mat
% 
% Hadi Hafizi, Nov. 2015

clear all

cd 091
JitBin = [2 5 10 15 25 30];

load Sim_bin1ms.mat
PairSim = S;
clear S
% S_r = PairSim(1:length(Rich_idx),1:length(Rich_idx));
% S_nr = PairSim(length(Rich_idx)+1:end,length(Rich_idx)+1:end);
% S_r_nr = PairSim(1:length(Rich_idx),length(Rich_idx)+1:end);

p_val = zeros(size(JitBin));
for ii=1:length(JitBin)
    load(['Sim_jit',num2str(JitBin(ii)),'bin_1ms.mat'])
%     S_shuf_r = mean(S_shuf(1:length(Rich_idx),1:length(Rich_idx),:),3);
%     S_shuf_nr = mean(S_shuf(length(Rich_idx)+1:end,length(Rich_idx)+1:end,:),3);
%     S_shuf_r_nr = mean(S_shuf(1:length(Rich_idx),length(Rich_idx)+1:end,:),3);
    p_val1 = zeros(size(S,3),1);
    parfor kk=1:100
    SimJit = S(:,:,kk);
    
    [h,p_val1(kk)] = kstest2(reshape(PairSim,size(PairSim,1)*size(PairSim,2),1),...
        reshape(SimJit,size(SimJit,1)*size(SimJit,2),1));
%     [h,p_val(ii)] = kstest2(reshape(PairSim,size(PairSim,1)*size(PairSim,2),1),...
%         reshape(SimJit,size(SimJit,1)*size(SimJit,2)*size(SimJit,3),1));
    end
    clear S
    p_val(ii) = mean(p_val1);
end


figure;
plot(JitBin,p_val,'LineWidth',2);hold on
plot(JitBin,0.05*ones(size(JitBin)),'LineWidth',2,'color','r');
title('KS-test for Similarity of Jittered Data','FontSize',16)
xlabel('Jittering value [ms]','FontSize',16);
ylabel('p-value','FontSize',16);
set(gca,'FontSize',16);
cd Figures
savefig('PairSimJitSigTest.fig')
cd ..


cd ..





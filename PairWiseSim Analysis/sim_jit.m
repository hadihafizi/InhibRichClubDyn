% How jittering effects similarity.
% 
% Required functions: sim_spktrn.m
% 
% Required files: asdf.mat
%                 PDF_1_16_30ms.mat
%                 wgts_1_16ms.mat
% 
% Sunny Nigam, Sept. 2015

clear all

load asdf.mat
load PDF_1_16_30ms.mat; load wgts_1_16ms.mat
W = PDF(:,:,45).*wgt;oute = sum(W,2);
[A2 B2] = sort(oute,'descend'); % sorted richness
duration = asdf_raw{end}(2);[A B] = size(W);
asdf_grp = ASDFSubsample(asdf_raw,B2); bin = 1;
asdf_bin = ASDFChangeBinning(asdf_grp,bin);
% S_act = sim_spktrn(asdf_bin);
jitter = [2 5 10 15 25 30];
for kk = 1:length(jitter);
    parfor ii = 1:100
           asdf_j = ASDFSpikeJittering(asdf_bin,jitter(kk),1);
           S(:,:,ii) = sim_spktrn(asdf_j);
   
    end
       str = ['save Sim_jit',num2str(jitter(kk)),'bin_',num2str(bin),'ms S;']; eval(str);
       clear S 
end
% str = ['save Sim_act_bin',num2str(bin),'ms S_act;']; eval(str);
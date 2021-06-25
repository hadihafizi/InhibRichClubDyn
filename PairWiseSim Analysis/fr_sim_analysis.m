
% idea similarity of connected and un-connected rich neurons
% re arrange the ASDF so that the nodes are ranked
% according to decreasing richness 
% 
% Required functions: sim_spktrn.m
%                     shuf_isi.m
% 
% Required files: asdf.mat
%                 PDF_1_16_30ms.mat
%                 wgts_1_16ms.mat
% 
% Sunny Nigam, Sept. 2015
% 
clear all
dataSet = '020';
load([dataSet,'\asdf.mat'])
% load([dataSet,'\PDF_NoSpur_thr45_',dataSet,'.mat']); 
load([dataSet,'\PDF_1_16_30ms.mat']); 
load([dataSet,'\wgts_1_16ms.mat'])
% W = PDF_cor.*wgt; 
W = PDF(:,:,45).*wgt;
oute = sum(weightmat,2);
inte = sum(weightmat,1);
tote = oute + inte';
[A2, B2] = sort(tote,'descend'); % sorted richness
% [A2, B2] = sort(oute,'descend'); % sorted richness
% duration = asdf_raw{end}(2);[A, B] = size(W);
asdf_grp = ASDFSubsample(asdf,B2);

bin = 1;
% bin = [1000 1500 2000 2500 3000 3500 4000 4500];
tic
for mm = 1:length(bin)
    asdf_bin = ASDFChangeBinning(asdf_grp,bin(mm));
    S = sim_spktrn(asdf_bin);
    parfor kk = 1:100
        [asdf_shuf] = shuf_isi(asdf_bin);
        %parsave(sprintf('asdf_shuf%d.mat', kk), asdf_shuf);
        S_shuf(:,:,kk) = sim_spktrn(asdf_shuf);
    end
%     Sshuf_mn(:,:,mm) = mean(S_shuf,3);
%     cor_sim(:,:,mm) = S - Sshuf_mn(:,:,mm);
    str = ['save Sim_bin',num2str(bin(mm)),'ms',' S S_shuf;']; eval(str); 
end
toc
%save('Sim_139.mat','S')

figure;
subplot(1,2,1);
imagesc(S+S');colorbar; colormap('jet');
subplot(1,2,2);
imagesc(mean(S_shuf,3)+transpose(mean(S_shuf,3)));colorbar; colormap('jet');
% subplot(1,3,3);
% cor_Sim = S-Sshuf_mn; cor_Sim(cor_Sim<0)=0;
%imagesc(cor_Sim);colorbar; colormap('jet');
figure; 
histogram(nonzeros(S),100,'normalization',...
   'probability'); hold on
histogram(nonzeros(S_shuf),100,'normalization',...
   'probability');
% save('1.mat','S','S_shuf')


    
    
    
    %SIM = sim+sim';
    %clearvars -except SIM oute_gr
    
%     [M Q]=community_louvain(SIM,1.1);
%     Q_max(mm) = Q;
%     comm_ind = unique(M);
%     for jj = 1:length(comm_ind)
%         com_mem{jj} = find(M == comm_ind(jj));
%         com_oute(jj) = mean(oute_gr(com_mem{jj}));
%         var_oute(jj) = std(oute_gr(com_mem{jj}));
%         temp = SIM(com_mem{jj},com_mem{jj});
%         temp1 = triu(temp);
%         sim_com(jj) = mean(nonzeros(temp1));
%         var_sim(jj) = std(nonzeros(temp1));
%     end
%     subplot(1,2,1)
%     errorbar(sim_com,var_sim);hold on
%     subplot(1,2,2)
%     errorbar(com_oute,var_oute);hold on
    %imagesc(SIM);colorbar; colormap('jet');axis square
    %addpath('C:\Users\Sunny\Desktop\Spiketrain_similarity\SpikeTrainCommunitiesToolBox-master')
%     tic
%     [grps,Qmax,grpscon,Qcon,ctr,maxQ,varargout] = allevsplitConTransitive(SIM);
%     maxmod(mm) = Qmax; max_cons(mm) = Qcon;
%     toc
%end
%figure;errorbar(bin,sim_val,var_sim)
%figure; plot(Q_max,'-o')
%plot(bin,maxmod,'-or'); hold on; plot(bin,max_cons,'-ob')
% % make synthetic spike trains with same ISIs and mean rate
% fr = ASDFGetfrate(asdf_raw).*1000;
% fr_r1 = fr(rich1); fr_r2 = fr(rich2);
% isi_r1 = diff(spk_r1); isi_r2 = diff(spk_r2);
% 
% tic
% for ii = 1:1000
%     syn_isi_r1 = datasample(isi_r1,length(isi_r1),'replace',false);
%     syn_isi_r2 = datasample(isi_r2,length(isi_r2),'replace',false);
%     % initialize first spike time
%     ind1 = 5;  % this is in ms
%     synspk_tr1 = ind1+cumsum(syn_isi_r1);
%     synspk_tr1 = horzcat(ind1,synspk_tr1);
% 
%     synspk_tr2 = ind1+cumsum(syn_isi_r2);
%     synspk_tr2 = horzcat(ind1,synspk_tr2);
% 
% 
% %generate asdf for synthetic spike trains
%     sub_r_syn = cell(4,1);
%     sub_r_syn{1} = synspk_tr1;
%     sub_r_syn{2} = synspk_tr2;
%     sub_r_syn{3} = 1;
%     sub_r_syn{4}(1,1) = 2; sub_r_syn{4}(1,2)= duration;
% 
% %Bin synthetic spike trains
%     bin_r_syn = ASDFChangeBinning(sub_r_syn,50);
% 
% % calulate similarity of above synthetic trains
%     inter_syn = intersect(bin_r_syn{1},bin_r_syn{2});
%     unio_syn = union(bin_r_syn{1},bin_r_syn{2});
%     sim_syn(ii) = length(inter_syn)/length(unio_syn);
%     clearvars -except sim_syn ii isi_r1 isi_r2 duration sim_a
% end
% toc
% histogram(sim_syn,50,'normalization','probability'); hold on
% plot([sim_a sim_a],[0 0.08],'-r','linewidth',2)





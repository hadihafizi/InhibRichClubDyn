function [str_aval] = Avalanche_detection(bs,asdf_raw)


% % detect avalanches from spike times of multiple neurons
% % Sunny Nigam 09/17/2015

% % 
% clear all; 
% cd(dataset)
% bs = 1; %bin size in [ms]
% addpath C:\Users\Hadi Hafizi\Box Sync\Avalanche_analysis
% load(['asdf_rest',dataset,'_TimeCorOneMSGap.mat']); %asdf_raw = asdf_shuf;
asdf_bin = ASDFChangeBinning(asdf_raw,bs);% clear asdf_raw
rast = ASDFToSparse(asdf_bin); clear asdf_bin
pop_fir = sum(rast);
%  find out where the population spike num~=0
ind_fr = find(pop_fir==0); clear pop_fir
                                                                                                       
% find time bins of continued activity
tic
kk = 1; 
for ii = 1:length(ind_fr)-1   
    if ind_fr(ii+1)-ind_fr(ii)>2  % leave out length 1 avalanche 
        ava = ind_fr(ii)+1:ind_fr(ii+1)-1;
        temp = rast(:,ava); clear ava % isolate those bins in the time raster
        [I J] = find(temp~=0); % find out which neurons spike and when
        sze(kk) = length(I); len(kk) = length(unique(J));
        % (relative times)
        str_aval{kk} = horzcat(I,J); clear temp I J % store this in a cell array
        kk = kk+1;
    end
    
end
toc
% save(['aval_list',num2str(dataset),'_bs',num2str(bs),'ms.mat'],'str_aval')
% cd ..
% each cell array in the above saved variable has 2 columns. The first
% column contains the neuron ID and the second column contains the time
% bin in which it spikes in the exact same sequence as the avalanche
% occured in the actual data. 

%plot length distribution and size distribution for all 
% [f1,f2] = hist(len,1:max(len)); [f3,f4] = hist(sze,1:max(sze));
% figure
% subplot(1,2,1);loglog(f2,f1,'linewidth',2); title('length distribution')
% set(gca,'linewidth',2,'fontsize',16)
% subplot(1,2,2);loglog(f4,f3,'linewidth',2); title('size distribution')
% set(gca,'linewidth',2,'fontsize',16)

% % make a list of rich and non rich nodes
% clear all
% load aval_list06-0_bs5ms.mat 
% load PDF_1_16_30ms.mat; load wgts_1_16ms.mat
% qwe = PDF(:,:,45).*wgt; oute = sum(qwe,2);
% [A2,B2]= sort(oute,'descend'); num_rich = ceil(0.2*length(A2));
% rich = B2(1:num_rich); 
% % nr1 = B2(num_rich+1:2*num_rich); % select the next set of nodes
% nr1 = B2(end-num_rich+1:end);
% % nr1 = B2(2*num_rich+1:3*num_rich);
% % (same number as the rich from the list)
% clear PDF mte_shuf mte fin_err cnd mte qwe max_ind wgt
% 
% % seperate out avalanches into 2 categories. 1st category contains
% % avalanches involving only rich nodes and the 2nd contains avalanches 
% % involving only non-rich nodes
% kk = 1; mm = 1;tic
% for ii = 1:max(size(str_aval))
%     mem = (unique(str_aval{ii}(:,1))); l1 = length(mem);
% %      if sum(ismember(mem,rich)) == l1 
%        if sum(ismember(mem,rich)) ~= l1
%          ava_rich{kk} = str_aval{ii};
%          len_aval_r(kk) = len(ii);sze_aval_r(kk) = sze(ii);
%          kk = kk +1;
% %      elseif sum(ismember(mem,nr1)) == l1 
% %        elseif sum(ismember(mem,nr1)) == l1 
% %          ava_nr1{mm} = str_aval{ii};
% %          len_aval_nr1(mm) = len(ii);sze_aval_nr1(mm) = sze(ii);
% %          mm = mm +1;
%      end
%        
% end
% toc
% 
% % length and size distribution of avalanches in rich sub-network
% [f1,f2] = hist(len_aval_r,1:max(len_aval_r)); 
% % [f3,f4] = hist(len_aval_nr1,1:max(len_aval_nr1)); 
% [f5,f6] = hist(sze_aval_r,1:max(sze_aval_r));
% % [f7,f8] = hist(sze_aval_nr1,1:max(sze_aval_nr1));
% figure
% subplot(1,2,1);loglog(f2,f1,'--or','linewidth',2); title('length distribution')
% % hold on; loglog(f4,f3,'--ob','linewidth',2);legend('rich','non-rich')
% set(gca,'linewidth',2,'fontsize',16)
% subplot(1,2,2);loglog(f6,f5,'--or','linewidth',2); title('size distribution')
% % hold on; loglog(f8,f7,'--ob','linewidth',2);legend('rich','non-rich')
% set(gca,'linewidth',2,'fontsize',16)
%  save('aval_listr_nr_both_bs5ms.mat','ava_rich','sze_aval_r','len_aval_r')%...
%      %'sze_aval_nr1',,'len_aval_nr1','ava_nr1')
% 
% % extract all avalanches of a particular length/size from each of these 
% % sub networks and calculate the similarity matrix for 
% % both. We want to check whether avalanches in the rich 
% % sub-network are more similar to each other than avalanches
% % in the non rich sub network 
% 
% clear all
% load('aval_listr_nrbs1ms.mat')
% len_chosen = [4:6]; sze_chosen = [3:5];
% % for ww = 1:length(len_chosen)
% %     tic
% %     indr = find(len_aval_r == len_chosen(ww));
% %     indnr = find(len_aval_nr1 == len_chosen(ww));
% %     
% % 
% %     for ii = 1:length(indr)
% %         list_r{ii} = ava_rich{indr(ii)}(:,1);
% %     end
% %     s1= max(size(list_r));
% % 
% %     for jj = 1:length(indnr)
% %         list_nr{jj} = ava_nr1{indnr(jj)}(:,1);
% %     end
% %     s2= max(size(list_nr));
% %     
% %     sim_r = sim_aval(list_r); sim_nr = sim_aval(list_nr);
% %     
% % %     figure;histogram(triu(sim_r),20); hold on; 
% % %     histogram(triu(sim_nr),20);
% %     S_R(ww) = sum(sum(sim_r))/(0.5*s1*(s1-1));
% %     S_NR(ww) = sum(sum(sim_nr))/(0.5*s2*(s2-1));
% % 
% %     clear list_r list_nr sim_r sim_nr indr indnr
% %   toc  
% % end
% 
% for ww = 1:length(sze_chosen)
%     tic
%     
%     indsr = find(sze_aval_r == sze_chosen(ww));
%     indsnr = find(sze_aval_nr1 == sze_chosen(ww));
%     
%     for ii = 1:length(indsr)
%         lists_r{ii} = ava_rich{indsr(ii)}(:,1);
%     end
%     s3= max(size(lists_r));
%     
%     for jj = 1:length(indsnr)
%         lists_nr{jj} = ava_nr1{indsnr(jj)}(:,1);
%     end
%     s4= max(size(lists_nr));
%     
% 
%     
%     sim_rs = sim_aval(lists_r); sim_nrs = sim_aval(lists_nr);
% 
%     S_Rs(ww) = sum(sum(sim_rs))/(0.5*s3*(s3-1));
%     S_NRs(ww) = sum(sum(sim_nrs))/(0.5*s4*(s4-1));
%     clear list_r list_nr sim_r sim_nr indr indnr
%   toc  
% end
% figure
% % subplot(1,2,1)
% % plot(len_chosen,S_R,'-or'); hold on; plot(len_chosen,S_NR,'-ob');
% % legend('rich','non-rich')
% % subplot(1,2,2)
% plot(sze_chosen,S_Rs,'-or','linewidth',2); hold on; 
% plot(sze_chosen,S_NRs,'-ob','linewidth',2);
% legend('rich','non-rich')
% set(gca,'linewidth',2,'fontsize',16)
% savefig('sim_r_nr_ava_345.fig')







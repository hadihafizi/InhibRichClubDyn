clear all

dataSet = {'02-0','03-0','03-1','04-0','04-1','05-0','05-1','06-0',...
    '06-1','07-0','07-1','08-0','08-1','09-0','09-1'};
suffix = '_cor_toteRC_Shuf';
% parpool(10);
for i = 1:length(dataSet)
%     asdf = load([dataSet{i},'\asdf.mat']);
    pdf = load([dataSet{i},'\PDF_NoSpur_thr45_',dataSet{i},'.mat']);
    wt = load([dataSet{i},'\wgts_1_16ms.mat']);
    % W = PDF_cor.*wgt;
    weightmat = (pdf.PDF).*(wt.wgt);
    oute = sum(weightmat,2);
    inte = sum(weightmat,1);
    tote = oute + inte';
    [A2, B2] = sort(tote,'descend'); % sorted richness
%     asdf_grp = ASDFSubsample(asdf.asdf_raw,B2);
    
    bin = 1;
    tic
%     asdf_bin = ASDFChangeBinning(asdf_grp,bin);
%     S = sim_spktrn(asdf_bin);
    S_shuf = zeros([size(weightmat),20]);
    parfor shufnum = 1:20
        asdf = load([dataSet{i},'\asdf_shuf',num2str(shufnum),'.mat']);
        asdf_grp = ASDFSubsample(asdf.x,B2);
%         [asdf_shuf] = shuf_isi(asdf_bin);
%         %parsave(sprintf('asdf_shuf%d.mat', kk), asdf_shuf);
        S_shuf(:,:,shufnum) = sim_spktrn(asdf_grp);
    end
    %     Sshuf_mn(:,:,mm) = mean(S_shuf,3);
    %     cor_sim(:,:,mm) = S - Sshuf_mn(:,:,mm);
    parsave( ['Sim_bin',num2str(bin),'ms_',dataSet{i},suffix,'.mat'],{'S_shuf'},{S_shuf});
%     parsave( ['Sim_bin',num2str(bin),'ms_',dataSet{i},suffix,'.mat'],{'S'},{S});
    toc
    %{
    [Rich, Rich_idx] = sort(tote,'descend');
    num_rich = ceil(0.2*length(Rich));
    nonRich_idx = Rich_idx(num_rich+1:end);
    Rich_idx = Rich_idx(1:num_rich);
    nr_oute = oute(nonRich_idx);
    [nr_vals, nr_idx] = sort(nr_oute,'descend');
    
    S_r = S(1:length(Rich_idx),1:length(Rich_idx));
    S_nr = S(length(Rich_idx)+1:end,length(Rich_idx)+1:end);
    S_r_nr = S(1:length(Rich_idx),length(Rich_idx)+1:end);
    
    %% Data Dist (Optimal Binning Scheme)
    
    [~, edges] = histcounts(S(triu(S,1)~=0));
    [rrcount, ~] = histcounts(S_r(triu(S_r,1)~=0), edges);
    [nrnrcount, ~] = histcounts(S_nr(triu(S_nr,1)~=0), edges);
    [nrrcount, ~] = histcounts(S_r_nr(triu(S_r_nr)~=0), edges);
    rrcount = rrcount./(sum(rrcount)*(edges(2) - edges(1)));
    nrnrcount = nrnrcount./(sum(nrnrcount)*(edges(2) - edges(1))); 
    nrrcount = nrrcount./(sum(nrrcount)*(edges(2) - edges(1))); 
    
    centers = movsum(edges, 2)/2;
    centers = centers(2:end);
    
    figure;
    plot(centers,rrcount,'r--.','MarkerSize',10); hold on
    plot(centers,nrnrcount,'b--.','MarkerSize',10); hold on
    plot(centers,nrrcount,'k--.','MarkerSize',10); hold on
    
    % legend('','R-R','','NR-NR','','R-NR');
    legend('R-R','NR-NR','R-NR');
    title('Jaccard Coefficients ProbDensity Distributions', 'FontSize',16)
    xlabel('Jaccard Coefficient');
    ylabel('Probability Density');
    set(gca, 'FontSize',16,'xscale','log');
%}
end








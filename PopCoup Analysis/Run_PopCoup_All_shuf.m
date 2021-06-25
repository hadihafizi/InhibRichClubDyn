% examine population coupling within and between sub-networks


clear all; warning off;
dataSet = {'05-0', '08-0'};
suffix = '_cor_toteRC';  %'rw' for rewired nets; '' for originals; '_cor' for
binSize = 1; %ms
% parpool(12);
for i = 1:length(dataSet)
    if strcmp(dataSet{i},'03-1')
        shufnums = 8:20;
    elseif strcmp(dataSet{i}, '05-0')
        shufnums = 4:9;
    elseif strcmp(dataSet{i}, '08-0')
        shufnums = 2:20;
    end
    for shufnum = shufnums
%         load([dataSet{i},'\asdf_rest',dataSet{i},'_TimeCorOneMSGap.mat'])
        load(['C:\Users\mhafizi\Box Sync\Research\Beggs Lab\Projects\Rich Club Dynamics Project\in vitro\Avalanche_analysis\',dataSet{i},'\asdf_shuf',num2str(shufnum),'.mat'])
        % asdf_raw = ASDFChooseTime(asdf_raw, 1.6e6, 2e6);
        load([dataSet{i},'\PDF_NoSpur_thr45_',dataSet{i},'.mat'])
        load([dataSet{i},'\wgts_1_16ms.mat'])
        W = PDF.*wgt;oute = sum(W,2);inte = sum(W,1);tote = inte + oute';
        [A, B] = sort(tote,'descend');
        num_rich = ceil(0.2*length(A));
        rich = B(1:num_rich);
        rest = B(num_rich+1:end);
        
        asdf_rich = ASDFSubsample(x,rich);
        asdf_rich = ASDFChangeBinning(asdf_rich,binSize);
        
        asdf_rest = ASDFSubsample(x,rest);
        asdf_rest= ASDFChangeBinning(asdf_rest,binSize);
        
        [A1 B1] = size(asdf_rich);[A2 B2] = size(asdf_rest);
        
        interval = 50; % # of time bins lag/lead on either side to calculate PopCoup
        Coup_r_r = zeros(length(rich),2*interval+1);
        Coup_r_nr = zeros(length(rich),2*interval+1);
        
        tic
        parfor kk = 1:length(rich);
            % Coup R --> R
            spktms = asdf_rich{kk};%clearvars -except asdf_raw spktms
            ind = setdiff([1:A1-2],kk);
            Coup_r_r(kk,:) = pop_coup(A1,asdf_rich,spktms,ind,interval);
            % Coup R --> NR
            ind = 1:length(asdf_rest)-2;
            Coup_r_nr(kk,:) = pop_coup(A1,asdf_rest,spktms,ind,interval);
        end
        toc
        
        Coup_nr_nr = zeros(length(rest),2*interval+1);
        Coup_nr_r = zeros(length(rest),2*interval+1);
        tic
        parfor kk = 1:length(rest);
            % Coup NR --> NR
            spktms = asdf_rest{kk};%clearvars -except asdf_raw spktms
            ind = setdiff([1:A2-2],kk);
            Coup_nr_nr(kk,:) = pop_coup(A2,asdf_rest,spktms,ind,interval);
            % Coup NR --> R
            ind = 1:length(asdf_rich)-2;
            Coup_nr_r(kk,:) = pop_coup(A1,asdf_rich,spktms,ind,interval);
        end
        toc
        
        save(['PopCoupAll_',dataSet{i},suffix,'_bs',num2str(binSize),'ms_shuf',num2str(shufnum),'.mat'],'Coup_r_r','Coup_r_nr','Coup_nr_r','Coup_nr_nr')
    end
end

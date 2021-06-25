
clear all

dataset = {'02-0','03-0','03-1','04-0','04-1','05-0','05-1',...
    '06-0','06-1','07-0','07-1','08-0','08-1','09-0','09-1'};
% dataset = {'139','151','152','163','165','168','174'};


for i = 1:length(dataset)
    load([dataset{i},'/PDF_NoSpur_thr45_',dataset{i},'.mat'])
    load([dataset{i},'/wgts_1_16ms.mat'])
    sig_te = PDF_cor.*wgt; %Adj_mat(Adj_mat>0) = 1;
    inte = sum(sig_te,1); oute = sum(sig_te,2);
    [inrichnss, inrichnssIdx] = sort(inte,'ascend'); 
    [outrichnss, outrichnssIdx] = sort(oute,'ascend');
%     num_rich = ceil(0.2*length(inrichnss));
%     inRichIdx = inrichnssIdx(1:num_rich);
%     innonRichIdx = inrichnssIdx(num_rich+1:end);
%     outRichIdx = outrichnssIdx(1:num_rich);
%     outnonRichIdx = outrichnssIdx(num_rich+1:end);
    
    inteCounts = zeros(length(unique(inte)),size(sig_te,1));
    outeCounts = zeros(length(unique(oute)),size(sig_te,2));
    for k = 1:size(wgt,1) % go thru individual neurons to take zeros out when counting in/out TE
        
        Adj_in = sig_te(:,k);
%         binEdges = (unique(inrichnss)); %binEdges(isinf(binEdges)) = binEdges(2) - 1;
        [inteCounts(:,k) , inteBins] = hist(Adj_in(Adj_in~=0), unique(inte));
        Adj_out = sig_te(k,:);
%         binEdges = (unique(outrichnss)); %binEdges(isinf(binEdges)) = binEdges(2) - 1;
        [outeCounts(:,k) , outeBins] = hist(Adj_out(Adj_out~=0), unique(oute));
        
    end
    [inX, inY] = meshgrid(log10(inteBins),1:length(inte));
    [outX, outY] = meshgrid(log10(outeBins),1:length(oute));
    
    inZ = (log10(inteCounts./repmat(sum(inteCounts),length(unique(inrichnss)),1)))';
%     inZ(isnan(inZ)) = -10; inZ(isinf(inZ)) = -10;
    outZ = (log10(outeCounts./repmat(sum(outeCounts),length(unique(outrichnss)),1)))';
%     outZ(isnan(outZ)) = -10; outZ(isinf(outZ)) = -10;
    
%     figure; surf(inX,inY,log10(inteCounts(:,inrichnssIdx))','EdgeColor','none'); view(2);
%     figure(1); subplot(3,5,i);
%     surf(inX,inY,inZ(inrichnssIdx,:),'EdgeColor','none'); view(2); hold on
%     colorbar
%     set(gca, 'XScale', 'log')
%     title(['Data set ',dataset{i}])
%     xlabel('log_{10}(inTE)')
%     ylabel('Neurons in Ascending Order of inTE')
    figure(3); subplot(2,4,i);
    imagesc(inZ(inrichnssIdx,:))
    colorbar
    title(['Data set ',dataset{i}])
    xlabel('inTE bins')
    ylabel('Neurons in Ascending Order of inTE')
%     figure; surf(outX,outY,log10(outeCounts(:,outrichnssIdx))','EdgeColor','none'); view(2);
%     figure(2); subplot(3,5,i);
%     surf(outX,outY,outZ(outrichnssIdx,:),'EdgeColor','none'); view(2); hold on
%     colorbar
%     set(gca, 'XScale', 'log')
%     title(['Data set ',dataset{i}])
%     xlabel('log_{10}(outTE)')
%     ylabel('Neurons in Ascending Order of outTE')
    figure(4); subplot(2,4,i);
    imagesc(outZ(outrichnssIdx,:))
    colorbar
    title(['Data set ',dataset{i}])
    xlabel('outTE bins')
    ylabel('Neurons in Ascending Order of outTE')

    clearvars -except i k dataset
end







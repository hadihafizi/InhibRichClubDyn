
% Compute the average avalanche length that each neuron initiates (if any).

clear all

% dataset = {'02-0','03-0','03-1','04-0','04-1','05-0','05-1','06-0','06-1','07-0','07-1',...
%     '08-0','08-1','09-0','09-1'};
% dataset = {'139', '151', '152', '163', '165', '168', '174'};
dataset = {'rand','vweakRC','weakRC','medRC1','medRC2','strongRC','rwNoRC'};
PFeig = [102,102,103,104,104,107,103];
binSize = 1; %ms
avalLen = cell(length(dataset),1);
avalLenMean = cell(length(dataset),1);
avalLenStd = cell(length(dataset),1);
avalSize = cell(length(dataset),1);
avalSizeMean = cell(length(dataset),1);
avalSizeStd = cell(length(dataset),1);
% parpool(12);
for i = 1:length(dataset)
    tic
%     asdfStruct = load([dataset{i},'/asdf.mat']);
%     asdfStruct = load([dataset{i},'/asdf_rest',dataset{i},'_TimeCorOneMSGap.mat']);
    asdfStruct = load(['Simk',num2str(PFeig(i)),'_Generic_',dataset{i},'vLowerSeeds.mat'], 'asdf');
%     avalStruct = load([dataset{i},'/aval_list',dataset{i},'_bs',num2str(binSize),'ms.mat']);
    str_aval = Avalanche_detection(binSize,asdfStruct.asdf);
    parsave(['aval_list',dataset{i},'_bs',num2str(binSize),'ms_vLowerSeeds.mat'], {'str_aval'}, {str_aval});
%     avalLen = cell(asdf{end}(1),1);
    
    for n = 1:asdfStruct.asdf{end}(1)
        nSamp = 0;
        for iaval = 1:length(str_aval)
            
            temp = str_aval{iaval};
            init =  temp(temp(:,2) == 1 | temp(:,2) == 2 | temp(:,2) == 3);
            neur = find(init == n);
            if neur
                avalLen{i}{n}(nSamp+1) = length(unique(temp(:,2)));
                avalSize{i}{n}(nSamp+1) = length(unique(temp(:,1)));
                nSamp = nSamp + 1; 
            end

            
        end
    end
    
    avalLenMean{i} = cellfun(@mean,avalLen{i});
    avalLenStd{i} = cellfun(@std,avalLen{i});
    avalSizeMean{i} = cellfun(@mean,avalSize{i});
    avalSizeStd{i} = cellfun(@std,avalSize{i});
    
    toc
end

parsave(['AvalInitiator_AvgSizeLen_allCBM.mat'], ...
    {'avalLenMean', 'avalLenStd', 'avalSizeMean', 'avalSizeStd'}, ...
    {avalLenMean, avalLenStd, avalSizeMean, avalSizeStd})


%% Plot Figures

for i = 1:length(dataset)
    
    pdfStruct = load([dataset{i},'/PDF_NoSpur_thr45_',dataset{i},'.mat']);
    wgtStruct = load([dataset{i},'/wgts_1_16ms.mat']);
    sig_te = wgtStruct.wgt.*pdfStruct.PDF_cor;
    inte = sum(sig_te,1)'; oute = sum(sig_te,2); tote = inte + oute;
    
    %%%%%%% plot figures for inTE
    
    [inrichnss, inrichnssIdx] = sort(tote,'ascend'); 
    [outrichnss, outrichnssIdx] = sort(oute,'ascend'); 
    
    [xData, yData] = prepareCurveData( log10(inrichnss), avalLenMean{i}(inrichnssIdx) );
    [fitresult, gof] = fit( xData, yData, 'poly1' );
    
    % Plot fit with data.
    figure;
    h = plot(fitresult,'k');hold on; set(h,'LineWidth',2)
    scatter(log10(inrichnss), avalLenMean{i}(inrichnssIdx),25,log10(outrichnss))
    colormap jet    
    xlabel('log_{10}(inTE)');
    ylabel('Avg Avalanche Length Initiated')
    
    %%%%%%%% plot same figures for outTE
    

    [xData, yData] = prepareCurveData( log10(outrichnss), avalLenMean{i}(outrichnssIdx) );
    [fitresult, gof] = fit( xData, yData, 'poly1' );
    
    % Plot fit with data.
    figure;
    h = plot(fitresult,'k');hold on; set(h,'LineWidth',2)
    scatter(log10(outrichnss), avalLenMean{i}(outrichnssIdx),25,log10(inrichnss))
    colormap jet    
    xlabel('log_{10}(outTE)');
    ylabel('Avg Avalanche Length Initiated')
    
end

%% Aggregate all data sets

inteAll = []; outeAll = [];
avalLenMeanAll = [];
avalSizeMeanAll = [];

for i = 7:7%1:length(dataset)
%     pdfStruct = load([dataset{i},'/PDF_NoSpur_thr45_',dataset{i},'.mat']);
%     wgtStruct = load([dataset{i},'/wgts_1_16ms.mat']);
    load(['Simk',num2str(PFeig(i)),'_Generic_',dataset{i},'vLowerSeeds.mat'], 'sig_te');
%     sig_te = wgtStruct.wgt.*pdfStruct.PDF;
    inteAll = [inteAll;sum(sig_te,1)']; outeAll = [outeAll;sum(sig_te,2)];
    
    avalLenMeanAll = [avalLenMeanAll,avalLenMean{i}];
    avalSizeMeanAll = [avalSizeMeanAll,avalSizeMean{i}];
end
toteAll = inteAll + outeAll;
[inrichnssAll, inrichnssIdxAll] = sort(inteAll,'ascend');
[outrichnssAll, outrichnssIdxAll] = sort(outeAll,'ascend');
[totrichnssAll, totrichnssIdxAll] = sort(toteAll,'ascend');

%%%%% totTE
figure;

% subplot(2,2,1);
[xData, yData] = prepareCurveData( (totrichnssAll), (avalLenMeanAll(totrichnssIdxAll)) );
[fitresult, gof] = fit( xData, yData, 'poly1' );
% Plot fit with data.
h = plot(fitresult,'k');hold on; set(h,'LineWidth',2)
scatter((totrichnssAll), (avalLenMeanAll(totrichnssIdxAll)),25)
% colormap jet
xlabel('Total TE');
ylabel('Avg Avalanche Length Initiated')
title({'Avg Avalanche Length Initiated- CBM Rewired No RC', ['R^2=',num2str(gof.rsquare),', polyCoef=',num2str(fitresult.p1)]})
% axis([1e-10 1e-2 0 100])
set(gca,'tickDir','out','color','none','fontname','arial','linewidth',1,'fontsize',12,'xscale','log'); box off;

figure;
% subplot(2,2,2);
[xData, yData] = prepareCurveData( (totrichnssAll), (avalSizeMeanAll(totrichnssIdxAll)) );
[fitresult, gof] = fit( xData, yData, 'poly1' );
% Plot fit with data.
h = plot(fitresult,'k');hold on; set(h,'LineWidth',2)
scatter((totrichnssAll), (avalSizeMeanAll(totrichnssIdxAll)),25)
% colormap jet
xlabel('Total TE');
ylabel('Avg Avalanche Length Initiated')
title({'Avg Avalanche Size Initiated- CBM Rewired No RC', ['R^2=',num2str(gof.rsquare),', polyCoef=',num2str(fitresult.p1)]})
% axis([1e-10 1e-2 0 100])
set(gca,'tickDir','out','color','none','fontname','arial','linewidth',1,'fontsize',12,'xscale','log'); box off;

%%
%%%%% outTE

[xData, yData] = prepareCurveData( log10(outrichnssAll), log10(avalLenMeanAll(outrichnssIdxAll)) );
[fitresult, gof] = fit( xData, yData, 'poly1' );

% Plot fit with data.
figure;
h = plot(fitresult,'k');hold on; set(h,'LineWidth',2)
scatter(log10(outrichnssAll), log10(avalLenMeanAll(outrichnssIdxAll)),25,log10(inrichnssAll))
colormap jet
xlabel('log_{10}(outTE)');
ylabel('log_{10}(Avg Avalanche Length Initiated)')
title('Avg Avalanche Length Initiated- All in vitro')









% Compute the average avalanche length that each neuron initiates (if any).

clear all

% dataset = {'02-0','03-0','03-1','04-0','04-1','05-0','05-1','06-0','06-1','07-0','07-1',...
%     '08-0','08-1','09-0','09-1'};
dataset = {'02-0'};
binSize = 1; %ms
avalLen = cell(length(dataset),1);
avalLenMeantemp = cell(length(dataset),1);
avalLenMean = cell(length(dataset),1);
% parpool(12);
for i = 1:length(dataset)
    asdfStruct = load([dataset{i},'/asdf.mat']);
    for shufnum = 1:20
    avalStruct = load([dataset{i},'/aval_list',dataset{i},'shuf',num2str(shufnum),'_bs',num2str(binSize),'ms.mat']);
%     avalLen = cell(asdf_raw{end}(1),1);
    
    for n = 1:asdfStruct.asdf_raw{end}(1)
        nSamp = 0;
        for iaval = 1:length(avalStruct.str_aval)
            
            temp = avalStruct.str_aval{iaval};
            init =  temp(temp(:,2) == 1 | temp(:,2) == 2 | temp(:,2) == 3);
            neur = find(init == n);
            if neur
                avalLen{i}{shufnum}{n}(nSamp+1) = length(unique(temp(:,2)));
                nSamp = nSamp + 1; 
            end

            
        end
    end
    
    avalLenMeantemp{i}{shufnum} = cellfun(@mean,avalLen{i}{shufnum});

    end
    avalLenMean{i} = cellfun(@mean,avalLenMeantemp{i});
end

%% Plot Figures
%{
for i = 1:length(dataset)
    
    pdfStruct = load([dataset{i},'/PDF_NoSpur_thr45_',dataset{i},'.mat']);
    wgtStruct = load([dataset{i},'/wgts_1_16ms.mat']);
    sig_te = wgtStruct.wgt.*pdfStruct.PDF;
    inte = sum(sig_te,1)'; oute = sum(sig_te,2);
    
    %%%%%%% plot figures for inTE
    
    [inrichnss, inrichnssIdx] = sort(inte,'ascend'); 
    
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
    
    [outrichnss, outrichnssIdx] = sort(oute,'ascend'); 

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
%}
%% Aggregate all data sets

inteAll = []; outeAll = [];
avalLenMeanAll = [];

for i = 1:length(dataset)
    pdfStruct = load([dataset{i},'/PDF_NoSpur_thr45_',dataset{i},'.mat']);
    wgtStruct = load([dataset{i},'/wgts_1_16ms.mat']);
    sig_te = wgtStruct.wgt.*pdfStruct.PDF;
    inteAll = [inteAll;sum(sig_te,1)']; outeAll = [outeAll;sum(sig_te,2)];
    
    avalLenMeanAll = [avalLenMeanAll,avalLenMean{i}];
end

[inrichnssAll, inrichnssIdxAll] = sort(inteAll,'ascend');
[outrichnssAll, outrichnssIdxAll] = sort(outeAll,'ascend');

%%%%% inTE

[xData, yData] = prepareCurveData( log10(inrichnssAll), log10(avalLenMeanAll(inrichnssIdxAll)) );
[fitresult, gof] = fit( xData, yData, 'poly1' );

% Plot fit with data.
figure;
h = plot(fitresult,'k');hold on; set(h,'LineWidth',2)
scatter(log10(inrichnssAll), log10(avalLenMeanAll(inrichnssIdxAll)),25,log10(outrichnssAll))
colormap jet
xlabel('log_{10}(inTE)');
ylabel('log_{10}(Avg Avalanche Length Initiated)')
title('Avg Avalanche Length Initiated- All in vitro')

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








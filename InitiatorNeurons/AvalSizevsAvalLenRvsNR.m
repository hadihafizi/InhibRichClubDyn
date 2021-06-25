clear all
% AvalSize vs AvalLen for low vs high TE neurons

dataset = {'02-0','03-0','03-1','04-0','04-1','05-0','05-1',...
    '06-0','06-1','07-0','07-1','08-0','08-1','09-0','09-1'};

load('AvgAvalLenInitiated.mat')
load('AvgandStdAvalSizeInitiated_Allinvitro.mat')
% inteAll = []; outeAll = [];
avalLenMeanRich = []; avalLenMeanNRich = [];
avalSizeMeanRich = []; avalSizeMeanNRich = [];
TopRichnssAll = []; BottRichnssAll = [];

for i = 1:length(dataset)
    load([dataset{i},'/PDF_NoSpur_thr45_',dataset{i},'.mat'])
    load([dataset{i},'/wgts_1_16ms.mat'])
    sig_te = PDF.*wgt; %Adj_mat(Adj_mat>0) = 1;
    inte = sum(sig_te,1); oute = sum(sig_te,2);
    [inrichnss, inrichnssIdx] = sort(inte,'ascend'); 
    [outrichnss, outrichnssIdx] = sort(oute,'ascend');
    num_rich = ceil(0.5*length(inrichnss));
    TopInIdx = inrichnssIdx(1:num_rich);
    BottInIdx = inrichnssIdx(end-num_rich+1:end);
    TopOutIdx = outrichnssIdx(1:num_rich);
    BottOutIdx = outrichnssIdx(end-num_rich+1:end);
    
    TopRichnssAll = [TopRichnssAll;outrichnss(1:num_rich)];
    BottRichnssAll = [BottRichnssAll;outrichnss(end-num_rich+1:end)];
    
    avalLenMeanRich = [avalLenMeanRich,avalLenMean{i}(TopOutIdx)];
    avalLenMeanNRich = [avalLenMeanNRich,avalLenMean{i}(BottOutIdx)];
    avalSizeMeanRich = [avalSizeMeanRich,avalSizeMean{i}(TopOutIdx)];
    avalSizeMeanNRich = [avalSizeMeanNRich,avalSizeMean{i}(BottOutIdx)];
    
    clearvars -except dataset i avalLenMeanRich avalLenMeanNRich avalSizeMeanRich avalSizeMeanNRich avalLenMean avalSizeMean TopRichnssAll BottRichnssAll
end

[xDataR, yDataR] = prepareCurveData( (avalLenMeanRich), (avalSizeMeanRich) );
[xDataNR, yDataNR] = prepareCurveData( (avalLenMeanNRich), (avalSizeMeanNRich) );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'off';
opts.Lower = [1 -Inf];
% opts.StartPoint = [0.629 1.58856397405898];
opts.Upper = [1 Inf];

[fitresultR, gofR] = fit( xDataR, yDataR, 'power1', opts );
[fitresultNR, gofNR] = fit( xDataNR, yDataNR, 'power1', opts );

% Plot fit with data.
figure;
h = plot(fitresultR,'r');hold on; set(h,'LineWidth',2)
h = plot(fitresultNR,'b');hold on; set(h,'LineWidth',2)

scatter((avalLenMeanRich), (avalSizeMeanRich),25, log10(TopRichnssAll))
scatter((avalLenMeanNRich), (avalSizeMeanNRich),25, log10(BottRichnssAll))
colormap jet
set(gca, 'xscale','log','yscale','log')
axis tight
l = legend(['Rich - T^{',num2str(fitresultR.b),'}'],...
    ['nonRich - T^{',num2str(fitresultNR.b),'}'],'Location','NorthWest');
set(l,'FontSize',12)
xlabel('Avg Avalanche Length Initiated');
ylabel('Avg Avalanche Size Initiated')

clear all

dataset = {'02-0','03-0','03-1','04-0','04-1','05-0','05-1','06-0','06-1','07-0','07-1',...
    '08-0','08-1','09-0','09-1'};
binSize = 1;
len = [2 5 8 12 15 20 25 30]; % desired length of avalanches to look at.
fig = 0; % whether to plot the figures

for i = 1:length(dataset)
    load([dataset{i},'/asdf.mat'])
    load([dataset{i},'/aval_list',dataset{i},'_bs',num2str(binSize),'ms.mat'])
    fr = cellfun(@numel, asdf_raw(1:end-2))/(asdf_raw{end}(2)/1000); % spikes/s
    
%     tic
%     [inits, nSamps] = AvalInits(str_aval,len);
%     toc
%     save([dataset{i},'/AvalInits_',dataset{i},'.mat'],'fr','inits','len','nSamps');
    load([dataset{i},'/AvalInits_',dataset{i},'.mat']);
    
    load([dataset{i},'/PDF_NoSpur_thr45_',dataset{i},'.mat'])
    load([dataset{i},'/wgts_1_16ms.mat']);
    sig_te = wgt.*PDF; oute = sum(sig_te,1);
    [richnss{i}, richnssIdx] = sort(oute,'descend'); num_rich = ceil(0.2*length(richnss{i}));
    if fig
%     figure(10*i+1);suptitle(['dataset: ',dataset{i}])
    H = figure(10*i+2);suptitle(['dataset: ',dataset{i}])
    set(H,'units','normalized','outerposition',[0 0 1 1])
    end
    for ilen = 1:length(len)
        [initCount, b] = histcounts(inits{ilen},asdf_raw{end}(1));
        initCount = initCount./fr'; % normalize by firing rate
        initsRSort{i}{ilen} = initCount(richnssIdx)/nSamps(ilen); % normalize by number of samples of each avalanche length
%         figure(10*i+1);
%         subplot(4,2,ilen);
%         semilogx(richnss, initsRSort, 'o', 'LineWidth', 2);
%         title(['Avalanche Length: ', num2str(len(ilen))]);
%         xlabel('TE Richness')
        
        %% Fit:
        if fig
        [xData, yData] = prepareCurveData( log10(richnss{i}), initsRSort{i}{ilen} );
        
        % Set up fittype and options.
        ft = fittype( 'poly1' );
        
        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft );
        
        % Plot fit with data.
        figure(10*i+2);subplot(4,2,ilen);
        h = plot( fitresult, xData, yData );
%         axis([-6.5 -2.5 0 .04])
        set(h(1),'LineWidth',2,'Color','k');set(h(2),'LineWidth',2,'Color','b')
        legend( h, 'Initiator-inRichness', 'Linear fit', 'Location', 'Best' );
        % Label axes
        xlabel('Log_{10}(TE inRichness)');
        ylabel({'#Avals norm. ','by nSamp and FR'})
        title(['Avalanche Length: ', num2str(len(ilen))]);
        end
        
    end
    clearvars -except i dataset binSize len richnss initsRSort fig
end

%% 
richnssAll = [];
for i = 1:length(dataset)
    richnssAll = [richnssAll,richnss{i}];
end
[richnssAll, richIdx] = sort(richnssAll,'ascend');

initsRSortAll = cell(length(len),1);
H = figure(10*i+5); %suptitle('All in vitro data sets')
set(H,'units','normalized','outerposition',[0 0 1 1])
for ilen = 1:length(len)
    for i = 1:length(dataset)
        initsRSortAll{ilen} = [initsRSortAll{ilen}, zscore(initsRSort{i}{ilen})];
    end
    initsRSortAll{ilen} = initsRSortAll{ilen}(richIdx);
    
    [xData, yData] = prepareCurveData( log10(richnssAll), initsRSortAll{ilen} );
    
    % Set up fittype and options.
    ft = fittype( 'poly1' );
    
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft );
    
    % Plot fit with data.
    figure(10*i+5);subplot(4,2,ilen);
    h = plot( fitresult, xData, yData );
    axis([-9 -2 -5 5])
    set(h(1),'LineWidth',2,'Color','k');set(h(2),'LineWidth',2,'Color','b')
    legend( h, 'Initiator-inRichness', 'Linear fit', 'Location', 'NorthWest' );
    % Label axes
    xlabel('Log_{10}(TE inRichness)');
    ylabel({'#Avals norm. ','by nSamp and FR'})
    title(['Avalanche Length: ', num2str(len(ilen))]);
        
end

%% PCA
%
i = 1;
% for i = 14:length(dataset)
    load([dataset{i},'/asdf.mat'])
    fr = cellfun(@numel, asdf_raw(1:end-2))/(asdf_raw{end}(2)/1000); % spikes/s
    load([dataset{i},'/AvalInits_',dataset{i},'.mat'],'inits','len');
    [initLongCount, ~] = histcounts(inits{8},asdf_raw{end}(1),'Normalization','pdf');
    [initLong, initLongIdx] = sort(initLongCount,'descend');
    
    
    load([dataset{i},'/PDF_NoSpur_thr45_',dataset{i},'.mat'])
    load([dataset{i},'/wgts_1_16ms.mat']);
    sig_te = wgt.*PDF; 
    inte = sum(sig_te,1)'; oute = sum(sig_te,2); 
    
    
% end



pcaVars = [fr,inte,oute];
pcaVars = bsxfun(@rdivide,pcaVars,mean(pcaVars));
[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(pcaVars);


pcagroups = zeros(asdf_raw{end}(1),1); pcagroups(initLongIdx(1:50)) = 1;
figure;
scatter3(SCORE(:,1),SCORE(:,2),SCORE(:,3),25,initLongCount,'filled')
title('PCA')
xlabel([num2str(COEFF(1,1)),'*fr + ',num2str(COEFF(2,1)),'*inTE + ',num2str(COEFF(3,1)),'*outTE'])
ylabel([num2str(COEFF(1,2)),'*fr + ',num2str(COEFF(2,2)),'*inTE + ',num2str(COEFF(3,2)),'*outTE'])
zlabel([num2str(COEFF(1,3)),'*fr + ',num2str(COEFF(2,3)),'*inTE + ',num2str(COEFF(3,3)),'*outTE'])

%}


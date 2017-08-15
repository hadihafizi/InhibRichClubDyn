clear all

dataSets = {'02-0','03-0','03-1','04-0','04-1','05-0','05-1','06-0',...
    '06-1','07-0','07-1','08-0','08-1','09-0','09-1'};

for i=1:length(dataSets)
    
    [FR, SpontFrac, InOutDeg, ISIstd] = InhibID(dataSets{i});
    
    figure(1);plot3k([ISIstd,SpontFrac,InOutDeg'],'Marker',{'.' 20});hold on %grid on
    ylabel('Fraction of Spontaneous Firing Rate')
    zlabel('Indeg - Outdeg')
    xlabel('std(ISI) [s]')

    figure(2);plot3k([FR,SpontFrac,InOutDeg'],'Marker',{'.' 20});hold on %grid on
    xlabel('Avg FR [spikes/s]')
    ylabel('Fraction of Spontaneous Firing Rate')
    zlabel('Indeg - Outdeg')

    figure(3);plot3k([FR,ISIstd,SpontFrac],'Marker',{'.' 20});hold on %grid on
    xlabel('Avg FR [spikes/s]')
    zlabel('Fraction of Spontaneous Firing Rate')
    ylabel('std(ISI) [s]')

    figure(4);plot3k([FR,ISIstd,InOutDeg'],'Marker',{'.' 20});hold on %grid on
    xlabel('Avg FR [spikes/s]')
    zlabel('Indeg - Outdeg')
    ylabel('std(ISI) [s]')
    
    clearvars -except dataSets i
end
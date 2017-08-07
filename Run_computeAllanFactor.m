
clear all

dataset = {'02-0','03-0','03-1','04-0','04-1','05-0','05-1','06-0','06-1','07-0','07-1',...
    '08-0','08-1','09-0','09-1'};

T = [50:50:1000, 1000:1000:1000000];
dt = 1;



AF = cell(length(dataset),1);%
% parpool(10);
for j = 1:1%length(dataset)
    load([dataset{j},'/asdf.mat'])
    AFtmp = zeros(length(T),asdf_raw{end}(1));
    for i = 1:asdf_raw{end}(1)
        spks = asdf_raw{i};
        Tstrt = asdf_raw{i}(1);
        Tend = asdf_raw{i}(end);
        AFtmp(:,i) = computeAllanFactor(spks, T, Tstrt, Tend, dt);
        
    end
    AF{j} = AFtmp;
end
save('AllanFactors_Allinvitro.mat','AF','T');

%% Plot Allan Factor for a Data Set

pdfStruct = load([dataset{1},'/PDF_NoSpur_thr45_',dataset{1},'.mat']);
wgtStruct = load([dataset{1},'/wgts_1_16ms.mat']);
sig_te = wgtStruct.wgt.*pdfStruct.PDF;
inte = sum(sig_te,1)'; oute = sum(sig_te,2);
[inrichnss, inrichnssIdx] = sort(inte,'ascend');
[outrichnss, outrichnssIdx] = sort(oute,'ascend');

[AFsort, AFsortIdx] = sort(sum(AF{j}),'ascend');

figure;
imagesc(AF{j}')
xlabel('Time window number')
ylabel('Neuron ID')

figure;
[Tmesh, Neurmesh] = meshgrid(T, 1:asdf_raw{end}(1));
surf(Tmesh, Neurmesh, log10(AF{j}(:,AFsortIdx))','EdgeColor','none');
view(2);
xlabel('Time Window [ms]')
ylabel('Neuron ID')
colormap jet

%% AF-based clusetering

j = 1;
AF{i}(isnan(AF{j})) = 0;
% AFsum = mean(AF{j});

nID = 1:asdf_raw{end}(1);
figure;
for i = 1:15
    AF_T = AF{j}(i,:)';
    idx = kmeans(AF_T,2);
    subplot(3,5,i);
    scatter(nID(idx == 1),(AF_T(idx == 1)),'r','LineWidth',2); hold on
    scatter(nID(idx == 2),(AF_T(idx == 2)),'b','LineWidth',2)
    xlabel('Neuron ID')
    ylabel('Allan Factor')
    title(['Tau = ',num2str(T(i)),' ms'])
end








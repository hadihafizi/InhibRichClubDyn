
clear all

dataset = {'02-0','03-0','03-1','04-0','04-1','05-0','05-1','06-0','06-1','07-0','07-1',...
    '08-0','08-1','09-0','09-1'};

load([dataset{1},'/asdf.mat'])
T = 5000:1000:1000000;
dt = 1;

pdfStruct = load([dataset{1},'/PDF_NoSpur_thr45_',dataset{1},'.mat']);
wgtStruct = load([dataset{1},'/wgts_1_16ms.mat']);
sig_te = wgtStruct.wgt.*pdfStruct.PDF;
inte = sum(sig_te,1)'; oute = sum(sig_te,2);
[inrichnss, inrichnssIdx] = sort(inte,'ascend');
[outrichnss, outrichnssIdx] = sort(oute,'ascend');

AF = zeros(length(T),asdf_raw{end}(1));
parpool(10);
parfor i = 1:asdf_raw{end}(1)
    spks = asdf_raw{i};
    Tstrt = asdf_raw{i}(1);
    Tend = asdf_raw{i}(end);
    AF(:,i) = computeAllanFactor(spks, T, Tstrt, Tend, dt);
    
end

[AFsort, AFsortIdx] = sort(sum(AF),'ascend');

figure;
imagesc(AF')
xlabel('Time window number')
ylabel('Neuron ID')

figure;
[Tmesh, Neurmesh] = meshgrid(T, 1:asdf_raw{end}(1));
surf(Tmesh, Neurmesh, log10(AF(:,AFsortIdx))','EdgeColor','none');
view(2);
xlabel('Time Window [ms]')
ylabel('Neuron ID')
colormap jet



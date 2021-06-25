% Having a list of avalanches, this code determines the richness profile 
% of avalanches with various lengths in a data set. It means the proportion 
% of active rich neurons are detected along the avalanches.
% 
% Hadi Hafizi and Sunny Nigam, Nov. 2015
% 
function [AvalLen,rchnss,err] = AvalRichProfileRev(len, str_aval, sig_te, varargin)

% if nargin > 3
%     inRich = varargin{1};
%     degMat = zeros(size(sig_te));
%     degMat(sig_te~=0) = 1;
%     if inRich
%         oute = sum(degMat,1);
%     else
%         oute = sum(degMat,2);
%     end
% end

oute = sum(sig_te,2);
inte = sum(sig_te,1);
tote = oute + inte';
% [A1,B1] = sort(oute,'descend'); 
[A1,B1] = sort(tote,'descend'); 
num_rich = ceil(0.2*length(A1));
% tmp = sum(sig_te,1);
% Exc = find(tmp>=0);
% rich = intersect(B1,Exc,'stable');
% rich = rich(1:num_rich);
rich = B1(1:num_rich);
% rich = B1(num_rich+1:end);

% average number of rich nodes involved in each time bin of 
% an avalanche of a particular length specified by len

avalnum = length(str_aval);
tic
for ilen = 1:length(len)
    str_cell{:} = [];
    qq = 1;
    for ii = 1:avalnum
        temp = str_aval{ii};
        L1 = length(unique(temp(:,2)));
        if L1 == len(ilen) || L1 == len(ilen)+1 || L1 == len(ilen)-1 || L1 == len(ilen)+2 || L1 == len(ilen)-2
%         if L1 == len(ilen) %|| L1 == len(ilen)+1 || L1 == len(ilen)-1 || L1 == len(ilen)+2 || L1 == len(ilen)-2
%             temp(:,2) = temp(:,2) - temp(1,2) + 1; % Rashid's way of saving avalanches is different than mine: he saves the actual timing of spikes but I start from 1 in all avalanches.
            for jj = 1:L1
                active = find(temp(:,2)==jj);
                active_rich = intersect(temp(active,1),rich);
                str(jj) = length(active_rich)/numel(active);
                clear nd_ri ind
            end
            str_cell{qq} = str;
%             str_cell(qq) = str;
            qq = qq+1;
            clear str
        end
        clear temp
    end
    % Collapsing AvalRichProfile vectors with different lengths to one
    % vector as an average of all samples
    [AvalLen{ilen},rchnss{ilen},err{ilen}] = AvgCollapse(str_cell,len(ilen));
    clear str_cell
    
end
toc
% 



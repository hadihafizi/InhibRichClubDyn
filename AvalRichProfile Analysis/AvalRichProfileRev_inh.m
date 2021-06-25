% Having a list of avalanches, this code determines the richness profile 
% of avalanches with various lengths in a data set. It means the proportion 
% of active rich neurons are detected along the avalanches.
% 
% Hadi Hafizi and Sunny Nigam, Nov. 2015
% 
function [AvalLen,rchnss,err, varargout] = AvalRichProfileRev_inh(len, str_aval, idx, varargin)

if nargin > 3
        idx1 = varargin{1};
        idx2 = varargin{2};
end        

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
                active_rich = intersect(temp(active,1),idx);
                str(jj) = length(active_rich)/numel(active);
                if nargin > 3
                    active_1 = intersect(temp(active,1),idx1);
                    str1(jj) = length(active_1)/numel(active);                    
                    active_2 = intersect(temp(active,1),idx2);
                    str2(jj) = length(active_2)/numel(active);                    
                end
                clear nd_ri ind
            end
            str_cell{qq} = str;
            if nargin > 3
                str_cell1(qq) = str1;
                str_cell2(qq) = str2;
                clear str1 str2
            end
            qq = qq+1;
            clear str
        end
        clear temp
    end
    % Collapsing AvalRichProfile vectors with different lengths to one
    % vector as an average of all samples
    [AvalLen{ilen},rchnss{ilen},err{ilen}] = AvgCollapse(str_cell,len(ilen));
    if nargin > 3
        [AvalLen{ilen},rchnss1{ilen},err1{ilen}] = AvgCollapse(str_cell1,len(ilen));
        [AvalLen{ilen},rchnss2{ilen},err2{ilen}] = AvgCollapse(str_cell2,len(ilen));
        varargout{1} = [AvalLen{ilen},rchnss1{ilen},err1{ilen}];
        varargout{2} = [AvalLen{ilen},rchnss2{ilen},err2{ilen}];
    end
    clear str_cell
    
end
toc
% 



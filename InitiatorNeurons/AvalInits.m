
function [inits, varargout] = AvalInits(str_aval, varargin)

% varargin: 
%           len: desired length of avalanches to inspect
%           nNeuron: number of neurons
% varargout:
%           nSamp: number of samples of each chosen avalanche lengths
%
% Hadi Hafizi, May 2017
%

if nargin > 1
    len = varargin{1};
    nSamp = zeros(length(len),1);
    inits = cell(length(len),1);
%     initNeurons = cell(length(len),1);
    for ilen = 1:length(len)
        for i = 1:length(str_aval)
            temp = str_aval{i};
            L1 = length(unique(temp(:,2)));
            if L1 == len(ilen) || L1 == len(ilen)+1 || L1 == len(ilen)-1 || L1 == len(ilen)+2 || L1 == len(ilen)-2
                inits{ilen} = [inits{ilen};temp(temp(:,2) == 1 | temp(:,2) == 2 | temp(:,2) == 3)];
                nSamp(ilen) = nSamp(ilen) + 1; % counting number of samples of each chosen Avalanche length
            end
        end
%         initNeurons{ilen} = unique(inits{ilen});
%         nNeuron = varargin{2};
%         [num(ilen), bins(ilen)] = hist(inits{ilen},0:nNeuron);
        if nargout > 1
            varargout{1} = nSamp;
        end
    end
else
    inits = [];
    for i = 1:length(str_aval)
        temp = str_aval{i};
        inits = [inits;temp(temp(:,2) == 1)];
    end
%     initNeurons = unique(inits);
%     [num, bins] = hist(inits,0:nNeuron);
end


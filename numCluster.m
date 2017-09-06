function numbClust = numCluster(data, varargin)

% Compute number of clusters in a data set using Sugar and James (2003)
% Jump method.
%
% Inputs: data:
%         maxK:
%         y:
%         plotjumps: Whether to plot figures (True/False)
%         printresults: Whether to print number of clusters found
%         (True/False)
% Output: numbClust:
% 
% Sugar, C. A., & James, G. M. (2003). Finding the number of clusters in a
% dataset: An information-theoretic approach. Journal of the American 
% Statistical Association, 98(463), 750-763.
% 
% Hadi Hafizi, Sept. 2017

% Define defaults
switch nargin
    case 1
        maxK = 10;
        y = [1, 1.5, 2];
        plotjumps = true;
        printresults = true;
    case 2
        maxK = varargin{1};
        y = [1, 1.5, 2];
        plotjumps = true;
        printresults = true;
    case 3
        maxK = varargin{1};
        y = varargin{2};
        plotjumps = true;
        printresults = true;
    case 4
        maxK = varargin{1};
        y = varargin{2};
        plotjumps = varargin{3};
        printresults = true;
    case 5
        maxK = varargin{1};
        y = varargin{2};
        plotjumps = varargin{3};
        printresults = varargin{4};
end
    
K = 1:maxK; % number of clusters to consider
n = size(data,1);
p = size(data,2);
sumd = zeros(length(K), 1);
for k = K
%   disp(['Clustering with K = ' num2str(k) ' ...']);
  [~, ~, sumdK] = kmeans(data, k, 'emptyaction', 'singleton');
  sumd(k) = sum(sumdK);
end
% compute distortion
dist = sumd/(n*p);

transDist = zeros(length(y),length(K)+1);
jumps = zeros(length(y),length(K));
numbClust = zeros(length(y),1);
for i = 1:length(y)
    % Compute the transformed distortion
    transDist(i,:) = [0, dist'.^(-y(i))];
    % Compute the jumps in transformed distortion
    jumps(i,:) = diff(transDist(i,:));
    % Compute the maximum jump
    numbClust(i) = K(jumps(i,:) == max(jumps(i,:)));
    % Plot distortion, transformed distortion and jumps
    if (plotjumps)
        figure(103);
        subplot(length(y),3,3*i-2);
        plot(K,dist,'LineWidth',2)
        xlabel('Number of Clusters')
        ylabel('Distortion')
        title(['Y = ',num2str(y(i))])
        subplot(length(y),3,3*i-1);
        plot([0,K],transDist(i,:),'LineWidth',2)
        xlabel('Number of Clusters')
        ylabel('Transformed Distortion')
        title(['Y = ',num2str(y(i))])
        subplot(length(y),3,3*i);
        plot(K,jumps(i,:),'LineWidth',2)
        xlabel('Number of Clusters')
        ylabel('Jumps')
        title(['Y = ',num2str(y(i))])
        % Plot line and point to indicate maximum jump
        
    end
    % Report maximum jump
    if (printresults)
        disp(['The maximum jump occurred at ',num2str(numbClust(i)), ' clusters with Y = ',num2str(y(i))]);
    end
end


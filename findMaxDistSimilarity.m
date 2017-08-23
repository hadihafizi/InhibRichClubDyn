function [ js_diverge ] = findMaxDistSimilarity( dat1, dat2, varargin )
% [ js_diverge, crossC ] = findMaxDistSimilarity( dat1, dat2 )
% Finds the Jensen-Shannon divergence between the probability distributions
% represented by two sets of data when one set is shifted so that the two
% possess maximal overlap. This function is meant to measure how different
% the shapes of two distributions are irrespective of their location along
% the x-axis. E.g. if the data in the first argument is drawn from a normal
% distribution with a mean of 0 and a standard deviation of 1 and the data
% in the second argument is drawn from a normal distribution with a mean of
% 2 and a standard deviation of 3, then this function will return the
% Jensen-Shannon divergence between the first set of data and the second
% set of data shifted along the x-axis by approximately -2, since that is
% where they maximally overlap. Another way to put this is that this
% function would return the same value as if the second data set where
% drawn from a normal disitrbution with a mean of 0 and a standard
% deviation of 2. The shift to creat maximum overlap is found by selecting
% the lag that maximized the cross-correlation between the two sets. Also
% returns the maximum cross correlation between the two binned datasets.
%
% [ js_diverge, crossC ] = findMaxDistSimilarity( dat1, dat2, 'Kernel', kernel )
%   Smooths the binned histograms of dat1 and dat2 before comparing them.
%   Can be used to get rid of noisy data or to make individual values in
%   bins not completely linearly independent when interpreted as a vector.
%
% [ js_diverge, crossC ] = findMaxDistSimilarity( dat1, dat2, 'Kernel', kernel, 'Bins', noBins)
%   Specify the number of bins to use when [dat1 dat2] is histogrammed.
%
% [ js_diverge, crossC ] = findMaxDistSimilarity( dat1, dat2, 'Kernel', kernel, 'Edges', edges)
%   Specify a vector of bin edges to use when binning [dat1 dat2]
%
% [ js_diverge, crossC ] = findMaxDistSimilarity( dat1, dat2, 'Kernel', kernel, 'BinMethod', binMethod)
%   Specify the binning method to use when bining [dat1 dat2]
%
% 'Bins', 'Edges', and 'BinMethod' are MUTUALLY EXCLUSIVE arguments. If
% more than one is used, the first will be used and the second ignored.
% This function will return an error if more than 2 optional arguments are
% used.

%==============================================================================
% Copyright (c) 2017, The Trustees of Indiana University
% All rights reserved.
%
% Authors:
% Z. Tosi (ztosi@indiana.edu),
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
%   1. Redistributions of source code must retain the above copyright notice,
%      this list of conditions and the following disclaimer.
%
%   2. Redistributions in binary form must reproduce the above copyright notice,
%      this list of conditions and the following disclaimer in the documentation
%      and/or other materials provided with the distribution.
%
%   3. Neither the name of Indiana University nor the names of its contributors
%      may be used to endorse or promote products derived from this software
%      without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%==============================================================================

narginchk(2, 9);

if ~isvector(dat1) || ~isvector(dat2)
    error('Data arguments must be vectors.');
end

binSpec = 0;
edSpec = 0;
binMethSpec = 0;

noBins = 0;
edges = [];
binMethod = '';
kernel = [];

ind = 1;
while ind <= length(varargin)
    switch varargin{ind}
        case 'Bins'
            if edSpec || binMethSpec
                warning('Bin edges or method have already been specified. Ignoring "Bin" argument.');
            else
                binSpec = 1;
                noBins = varargin{ind+1};
                if ~isscalar(noBins)
                    error('Bin # argument must be scalar!');
                end
            end
            ind = ind + 2;
        case 'Edges'
            if binSpec || binMethSpec
                warning('Number of bins or bin method have already been specified. Ignoring "Edges" argument.');
            else
                edSpec = 1;
                edges = varargin{ind+1};
                if ~isvector(edges)
                    error('Bin edges must be a vector!');
                end
            end
            ind = ind + 2;
        case 'BinMethod'
            if binSpec || edSpec
                warning('Number of bins or bin edges have already been specified. Ignoring "BinMethod" argument');
            else
                binMethSpec = 1;
                binMethod = varargin{ind+1};
            end
            ind = ind + 2;
        case 'Kernel'
            kernel = varargin{ind+1};
            ind = ind + 2;
        case 'Cutoff'
            dat1 = dat1(dat1<varargin{ind+1});
            dat2 = dat2(dat2<varargin{ind+1});
            ind = ind + 2;
        case 'Normalize'   
            dat1 = (dat1-mean(dat1))./std(dat1);
            dat2 = (dat2-mean(dat2))./std(dat2);
            ind=ind+1;
        otherwise
            error('Unrecognized input option');
    end
end

if isrow(dat1)
    dat1=dat1';
end
if isrow(dat2)
    dat2=dat2';
end

% Create bin-edges based on provided args...
if binMethSpec
    [~, edges] = histcounts([dat1; dat2], 'BinMethod', binMethod);
elseif binSpec
    [~, edges] = histcounts([dat1; dat2], noBins);
elseif edSpec
    % Do nothing because we already have edges...
else
    [~, edges] = histcounts([dat1; dat2]);
end

d1 = histcounts(dat1, edges, 'Normalization', 'probability');
d2 = histcounts(dat2, edges, 'Normalization', 'probability');

if ~isempty(kernel)
    % Smooth distributions with kernel
    d1 = conv(d1, kernel, 'same');
    d2 = conv(d2, kernel, 'same');
    % renormalize
    d1 = d1 ./ sum(d1);
    d2 = d2 ./ sum(d2);
end

d1=d1';
d2=d2';

% Find max cross correlation and lag where it happens...
[~, mxLag] = max(conv(d1, flipud(d2), 'full'));

% how much to shift d1
shift = mxLag - length(d1);

% Padded vectors for comparing shifted dist to non shifted dist,
% represent probability distributions p and q to be compared
% respectivly
p_pr = zeros(3*length(d1)-1, 1);
q_pr = zeros(3*length(d1)-1, 1);

% Indices in padded vector for d1 (unshifted) and d2 (shfited)
% respectivly
inds1 = (1:length(d1)) + length(d1)-1;
inds2 = inds1 + shift;

% Populate padded vectors.
p_pr(inds1) = d1;
q_pr(inds2) = d2;

%figure; plot(p_pr); hold on; plot(q_pr); hold off;

% Compute Jensen-Shannon Divergence in base-2
m_pr = (p_pr + q_pr)./2;
pnz = p_pr ~= 0;
qnz = q_pr ~= 0;
p_pr = 0.5 * sum(p_pr(pnz) .* log2(p_pr(pnz) ./ m_pr(pnz)));
q_pr = 0.5 * sum(q_pr(qnz) .* log2(q_pr(qnz) ./ m_pr(qnz)));
js_diverge = p_pr + q_pr;
end


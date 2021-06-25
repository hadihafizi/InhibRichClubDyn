function asdf_jitter = ASDFSpikeJittering(asdf, width, mode)
% asdf_jitter = ASDFSpikeJittering(asdf, width, mode)
%
%    asdf - {nNeu+2, 1} ASDF
%    width - (1,1) jittering width
%    mode - (1,1) 0 - normal dist 1 - uniform dist
%
% Returns:
%   asdf_jitter - {nNeu+2, 1} ASDF with jittered time
%
%
% Description :
%   Do spike jittering based on John's code
%
% Example :
%    
%
% Author   : Shinya Ito
%            Indiana University
%
% Last modified on 1/24/2011

if nargin < 3
	mode = 1;
end

raster = ASDFToSparse(asdf);
if mode == 0
	jittered = jitter2(raster, width);
else
	jittered = jitter2_u(raster, width);
end
asdf_jitter = SparseToASDF(jittered, asdf{end-1});

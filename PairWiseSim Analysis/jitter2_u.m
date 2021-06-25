function [JRASTER] = jitter2_u(RASTER, width)
%takes a RASTER (which may contain amplitude or binary values) and jitters the time
%of each event by a random number of bins drawn from a uniform distribution with width width.
%Note: some loss of events may take place here, since two events may be jittered into the same space, but this will be rare (I think)
%The jittered raster is returned as JRASTER.
% 01/21/02 JMB

%caution: may run out of memory on MEAN, so best to run this on MEAT
%if it must be run on a computer with small memory, you can divide RASTER into equally sized chunks, jitter them, then concatenate.

[r c] = size(RASTER);
[e f a] = find(RASTER);
J = round(width * rand(length(f), 1) - width/2);
Jf = f + J;
Jf(find(Jf <= 0)) = 1;
Jf(find(Jf > c)) = c;
JRASTER = sparse(e, Jf, a, r, c);

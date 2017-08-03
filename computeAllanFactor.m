function [AF] = computeAllanFactor(spks, T, Tstrt, Tend, dt)

% Caulculates Allan Factor based on "Teich et al 1997"
%
% Inputs: spks: vector of time stamps of spike times
%         T: vector of time-windows over which AF is calculated (in bins). 
%         (Tstrt - Tend)/T(i) should be larger than 2.
%         Tstrt: start time of data (spks) (scalar, in bins)
%         Tend: end time of data (scalar, in bins)
%         dt: bin size (in ms)
% 
% Output: AF: Allan Factor
%
% Version 0.1 (non-overlapping time windows)
%
% Hadi Hafizi, Aug. 2017
%
%

spks = spks - spks(1); % to remove any bias and make time stamps start from 0

Avar = zeros( length(T), 1 ); AF = zeros( length(T), 1 );

for i = 1:length(T)
    Avartmp = zeros( floor( ( Tend - Tstrt )/T(i) ) - 1, 1 ); 
    eventnum = zeros( floor( ( Tend - Tstrt )/T(i) ) - 1, 1 );
    for j = 1:floor( ( Tend - Tstrt )/T(i) ) - 1
        Avartmp(j) = ( length( spks( ( spks >= j*T(i)) & ( spks < ( j + 1)*T(i) ))) - ...
            length( spks( (spks >= (j - 1)*T(i)) & (spks < j*T(i) )))).^2;
        eventnum(j) = length( spks( (spks >= (j - 1)*T(i)) & (spks < j*T(i) )));
    end
    Avar(i) = mean(Avartmp);
    AF(i) = Avar(i) ./ ( 2*mean( eventnum ) );
    
end


figure;
plot( T*dt, AF, 'LineWidth', 2 )
title( 'Allan Factor' )
xlabel( 'Time Intervals [ms]' )
ylabel ( 'Allan Factor' )




function [TIMERASTERswapped] = spikeSwapper4(TIMERASTER)


% I. Identify spike times for each active channel
[r c] = size(TIMERASTER);
SumRows = sum(TIMERASTER,2);
ActiveChans = find(SumRows > 0);
TIMES = zeros(length(ActiveChans), max(SumRows));
for i=1:length(ActiveChans)
    [chans, times, s] = find(TIMERASTER(ActiveChans(i),:));
    TIMES(i,1:length(times)) = times;
end


% II. Arrange times so that they can be permuted; also account for
% different numbers of spikes in each train before permuting.
[Y I] = sort(SumRows(ActiveChans));
sortedTIMES = TIMES(I, :);
interval{1} = 1:Y(1);
for i=1:length(Y)-1
    interval{i+1} = Y(i)+1:Y(i+1);  
end



% III. Permute times (this will accomplish spike swapping)
swappedTIMES = sortedTIMES;
for i=1:length(ActiveChans)-1
    if ~isempty(interval{i})
        for j=interval{i}(1):interval{i}(length(interval{i}))
            chans = i:length(ActiveChans);
            swappedTIMES(chans, j) = sortedTIMES(chans(randperm(length(chans))), j);
        end
    end
end


% IV. Insert swapped spikes into TIMERASTERswapped
ColInd = zeros(sum(SumRows),1);
RowInd = zeros(sum(SumRows),1);
startPt = 1;
stopPt = SumRows(ActiveChans(I(1)));
for i=1:length(ActiveChans)
    ColInd(startPt:stopPt) = ActiveChans(I(i)) * ones(length(startPt:stopPt), 1);
    RowInd(startPt:stopPt) = swappedTIMES(i,1:SumRows(ActiveChans(I(i))));
    if i < length(ActiveChans)
        startPt = stopPt + 1;
        stopPt = (startPt - 1) + SumRows(ActiveChans(I(i+1)));
    end
end
TIMERASTERswapped = sparse(ColInd, RowInd, ones(length(sum(SumRows)), 1), r, c);
TIMERASTERswapped(find(TIMERASTERswapped > 1)) = 1;




% V. Correct for lost spikes
% Subtract swapped spike trains from original spike trains.
% All the positive values will be locations where a spike could be added to
% the swapped spike train. Insert these "real" spikes back into
% TIMERASTERswapped to make up for the lost spikes.
deltaTIMERASTER = TIMERASTER - TIMERASTERswapped;
RowInds = [];
ColInds = [];
for i=1:length(ActiveChans)
    numberMissing = full(sum(TIMERASTER(ActiveChans(i),:)) - sum(TIMERASTERswapped(ActiveChans(i),:)));
    [Z K] = find(deltaTIMERASTER(ActiveChans(i), :) > 0);
    Inds = randperm(length(K));
    %TIMERASTERswapped(ActiveChans(i), K(Inds(1:numberMissing))) = 1;
    RowInds = [RowInds; (ActiveChans(i) * ones(numberMissing, 1))];
    ColInds = [ColInds; K(Inds(1:numberMissing))'];
end
TIMERASTERswappedMissingSpikes = sparse(RowInds, ColInds, 1, r, c);
TIMERASTERswapped = TIMERASTERswapped + TIMERASTERswappedMissingSpikes;



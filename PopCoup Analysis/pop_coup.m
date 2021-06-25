function [tot_count] = pop_coup(A,asdf_bin,spktms,ind,interval)
% calculates how much the rest of the network spikes
% when a particular neuron fires in a -interval [ms] to +interval [ms]
% window. This is termed as population coupling.

    if length(spktms)>1
        parfor ii = 1:length(spktms);
            tmstr = spktms(ii)-interval; tmend = spktms(ii)+interval;
            newasdf = ASDFChooseTime(asdf_bin,tmstr,tmend);
            modasdf = ASDFSubsample(newasdf,ind);
            modasdf{end}=[]; modasdf{end-1}=[];
            bn_edge = 1:(2*interval)+1;
            [a, b]=hist(cell2mat(cellfun(@(x)x(:),modasdf(:),'un',0)),bn_edge);
            count(ii,:) = a; %clear a b       
        end
        tot_count  = sum(count)/length(spktms); %normalizes the count by 
        %the number of spikes in the reference neuron
        %plot([-500:500],tot_count); hold on
    elseif length(spktms) == 1
        for ii = 1:length(spktms);
            tmstr = spktms(ii)-interval; tmend = spktms(ii)+interval;
            newasdf = ASDFChooseTime(asdf_bin,tmstr,tmend);
            modasdf = ASDFSubsample(newasdf,ind);
            modasdf{end}=[]; modasdf{end-1}=[];
            bn_edge = 1:(2*interval)+1;
            [a, b]=hist(cell2mat(cellfun(@(x)x(:),modasdf(:),'un',0)),bn_edge);
            count(ii,:) = a; %clear a b       
        end
        tot_count = count;
    else
        tot_count = zeros(1,2*interval+1);
    end
    
   
end


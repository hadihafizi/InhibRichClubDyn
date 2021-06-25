% This function gets a cell array of proportion of active rich neurons in
% avalanches of a particular duration +/- 2 time steps and the
% corresponding duration. The outputs are "Un" that is a vector with the 
% points at which the "av" is calculated (x-axis). "av" is the ensemble 
% average proportion of active rich neurons in which vectors with different
% lengths are collapsed together. "err" is the standard error of the
% measure.
% 
% Hadi Hafizi, Jan 2016


function [Un,av,err] = AvgCollapse(AvalRichProportionCell,Length)

x1 = 1:(Length - 2);
x2 = 1:(Length - 1);
x3 = 1:Length;
x4 = 1:(Length + 1);
x5 = 1:(Length + 2);

x1scaled = (x1 - 1)/(length(x1) - 1);
x2scaled = (x2 - 1)/(length(x2) - 1);
x3scaled = (x3 - 1)/(length(x3) - 1);
x4scaled = (x4 - 1)/(length(x4) - 1);
x5scaled = (x5 - 1)/(length(x5) - 1);
% figure; hold on; plot(x1scaled,s1); plot(x2scaled,s2); plot(x3scaled,s3)
Un = union(union(union(x1scaled,union(x2scaled,x3scaled)),x4scaled),x5scaled);
% Sample
% s1 = -0.25*(x1 - 1).*(x1 - 10) + rand([1,10]);
% s2 = -0.25*(x2 - 1).*(x2 - 14) + rand([1,14]);
% s3 = -0.25*(x3 - 1).*(x3 - 18) + rand([1,18]);
% figure; hold on; plot(x1,s1); plot(x2,s2); plot(x3,s3)
s1int = []; s2int = []; s3int = []; s4int = []; s5int = [];
kk = 1; jj = 1; ww = 1; qq =1; pp = 1;
for ii = 1: length(AvalRichProportionCell)
    L1 = length(AvalRichProportionCell{ii});
    switch L1
        case Length - 2
            s1(kk,:) = AvalRichProportionCell{ii};
            s1int(kk,:) = interp1q(x1scaled',s1(kk,:)',Un')'; % My version of Matlab wants inputs to interp1q as columns
            kk = kk +1;
        case Length - 1
            s2(jj,:) = AvalRichProportionCell{ii};
            s2int(jj,:) = interp1q(x2scaled',s2(jj,:)',Un')';
            jj = jj +1;
        case Length
            s3(ww,:) = AvalRichProportionCell{ii};
            s3int(ww,:) = interp1q(x3scaled',s3(ww,:)',Un')';
            ww = ww +1;
        case Length + 1
            s4(qq,:) = AvalRichProportionCell{ii};
            s4int(qq,:) = interp1q(x4scaled',s4(qq,:)',Un')';
            qq = qq + 1;
        case Length + 2
            s5(pp,:) = AvalRichProportionCell{ii};
            s5int(pp,:) = interp1q(x5scaled',s5(pp,:)',Un')';
            pp = pp + 1;
    end
end
Un = Un*Length;
av = mean([s1int;s2int;s3int;s4int;s5int],1);
err = std([s1int;s2int;s3int;s4int;s5int],0,1)/...
    sqrt(length([s1int;s2int;s3int;s4int;s5int]));
% plot(Un,av,'r')
% figure;
% shadedErrorBar(Un,av,err,'r');hold on

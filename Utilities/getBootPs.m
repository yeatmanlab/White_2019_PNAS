%function pvals = getBootPs(bootDist, baseline)
%Little function to get pvalues from a bootstrapped distribution
%compared against some baseline value
%
%
%Inputs:  
%- bootDist: a BxC matrix of bootstrapped distributions, with B individual bootstrap samples in rows
%and C conditions in columns 
%- baseline: baseline level to compare against 
%
%Output: 
%- pvals, a 1 x C vector of p-values for a two-tailed test 
% 
% by Alex L. White at the University of Washington 

function pvals = getBootPs(bootDist, baseline)

B = size(bootDist,1); 
C = size(bootDist,2);

%get a large set of quantiles
qlevs=(1/(B+1)):(0.0001):(1-1/(B+1));
quants=quantile(bootDist,qlevs')';

%for each condition, find the quantile the baseline is closest to
pvals=NaN(1,C);
for ci=1:C
    qbasei=find(abs(quants(ci,:)-baseline)==min(abs(quants(ci,:)-baseline)));
    if ~isempty(qbasei)
        pvals(ci)=mean(qlevs(qbasei)); %mean in case there are 2 that perfectly straddle the baseline
    end
end
pvals(pvals>.5)=1-pvals(pvals>.5); %because some might be below baseline, so pval>.5;
pvals=2*pvals; %because this is a two-tailed test
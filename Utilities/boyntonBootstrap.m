function [CI,sampleStat,bootstrapStat] = boyntonBootstrap(myStatistic,x,nReps,CIrange,BCFlag)
% [CI,sampleStat] = bootstrap(myStatistic,x,[nReps],[CIrange],[BCFlag])
%
% Calculates a confidence interval on the statistic (function handle
% 'myStatistic' using a nonparametric bootstrap and the bias-corrected
% accellerated' algorithm.
%
% Inputs:
%    myStatistic:       a handle to a function that takes in a vector and
%                       returns a scalar (like 'mean' or 'std')
%
%   x                   sample vector for generating the statistic and its
%                       confidence intervals
%
%   nReps               number of sample-with-replacement iterations
%                       (default 2000)
%  
%   CIrange             confidence interval range (default 68.27)
%
%   BCFlag              boolean flag for whether or not to do the
%                       bias-correction (default 1)
%
% Outputs:     
%   CI:                 confidence interval
%   sampleStat          statictic evaluated on x  'myStatistic(x)'
%   bootstrapStat       vector of bootstrapped values of statistic   (1xnReps)

% 4/10/09 Written by G.M. Boynton at the University of Washington
%         based on  Efron and Tibshirani's "Introductcion to the Bootsrap" 
%        (Chapman & Hall/CRC, 1993, pages 178-201)
% 04/12/18 Edited by ALW to return bootstrapStat

%%
% Deal with defaults

if ~exist('BCFlag','var')
    BCFlag = 1;
end

if ~exist('CIrange','var')
    CIrange = 68.27;  %corresponds to +/- 1 s.d. for a normal distribution
end

if ~exist('nReps','var')
    nReps=2000;
end

%%
% Run the statistic on the sample 'x'

sampleStat = myStatistic(x);

%% 
% Run the statistic on resampled versions of x
id = ceil(rand(length(x),nReps)*length(x)); %matrix of random indices into x
resampledX = x(id);      %resampled versions of 'x

bootstrapStat = zeros(1,nReps);
for i=1:nReps
    bootstrapStat(i) = myStatistic(resampledX(:,i));
end


%% Calculate the bias-corrected and accelerated parameters z0 and a

if BCFlag
    z0 = norminv(sum(bootstrapStat<sampleStat)/nReps);

    thetai = zeros(1,length(x));
    %calculate the statistic holding one member of x out each time
    for i=1:length(x)
        id = [1:(i-1),(i+1):length(x)];
        thetai(i) = myStatistic(x(id));
    end
    %do something related to skewness.
    a = sum( (mean(thetai)-thetai).^3)/(6*(sum( (mean(thetai)-thetai).^2).^(3/2)));
else
    z0 =0;
    a=0;
end

%% Calculate the 'bias-corrected and accelerated' percentiles using z0 and a

zLo = norminv((1-CIrange/100)/2);
zHi = norminv((1+CIrange/100)/2);

zClo = z0 + (z0+zLo)/(1-a*(z0+zLo));
bcaLo = normcdf(zClo,0,1);

zChi = z0 + (z0+zHi)/(1-a*(z0+zHi));
bcaHi = normcdf(zChi,0,1);


CI(1) = prctile(bootstrapStat,100*bcaLo);
CI(2) = prctile(bootstrapStat,100*bcaHi);

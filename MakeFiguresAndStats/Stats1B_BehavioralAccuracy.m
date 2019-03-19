%% function Stats1B_BehavioralAccuracy() 
% Print stats about behavioral performance in White, Palmer, Boynton & Yeatman, PNAS 2019. 
% 
% by Alex L. White
% University of Washington, 2019

function Stats1B_BehavioralAccuracy()

%% set paths
%add whole analysis code directory to the path
analysisDir = fileparts(fileparts(which(mfilename)));
addpath(genpath(analysisDir));

p = getPaths();

%% load data
resFile = fullfile(p.data,'AllSubjBehavior.mat');
load(resFile);

%% stats file
statsFile = fullfile(p.stats,'Stats1B_AccuracyEffects.txt');
statsF = fopen(statsFile,'w');

fprintf(statsF,'Effects of cue and side on behavioral accuracy (p(correct)) \n\n');
%% extract data

conds.cueCond = [2 3]; %divided, focal
conds.targSide = 2:3; %all, left, right
conds.congruent = 1; %not divided by congruency
conds.targPres = 1; %not divided by target presence/absence

cueLabels = {'Distributed','Focal'};
pcs = squeeze(allR.PCs(conds.cueCond, conds.targSide, conds.congruent, conds.targPres,:));

%% cue effects
cueEffects = squeeze(diff(pcs,1,1));

%collapse across side
cueEffects = squeeze(mean(cueEffects,1));

[~,tp,~,tst] = ttest(cueEffects);
cueEffectCI = boyntonBootstrap(@mean,cueEffects,5000,95);

fprintf(statsF,'Focal - distributed cue effect (averaging across sides):\n');
fprintf(statsF,'\tmean = %.4f, SEM = %.4f, 95%% bootstrapped CI = [%.4f %.4f]\n\tt(%i) = %.4f, p=%.8f\n\n', mean(cueEffects), standardError(cueEffects), cueEffectCI(1), cueEffectCI(2), tst.df, tst.tstat, tp);

%% side effects
sideEffects = squeeze(diff(pcs,1,2));
for ci = 1:2
    fprintf(statsF,'\nRight - left side effect, in %s condition:\n', cueLabels{ci});
    meanEff = mean(sideEffects(ci,:));
    semEff = standardError(sideEffects(ci,:));
    CIEff = boyntonBootstrap(@mean, sideEffects(ci,:), 5000, 95);
    [~,tp,~,tst] = ttest(sideEffects(ci,:));
    fprintf(statsF,'\tmean = %.4f, SEM = %.4f, 95%% bootstrapped CI = [%.4f %.4f]\n\tt(%i) = %.4f, p=%.8f, ', meanEff, semEff, CIEff(1), CIEff(2), tst.df, tst.tstat, tp);
end

%% Serial model parameter estimates from each observer's AOC

fprintf(statsF,'\n\n\n------------------------\nSERIAL MODEL ESTIMATES\n');

p1s = allR.pProcessTask1First_PCs;
fprintf(statsF,'\nMean p(left side first) = %.2f +/- %.2f\n\t', mean(p1s), standardError(p1s));
[~,pval,~,sts] = ttest(p1s,0.5);
CI = boyntonBootstrap(@mean, p1s, 5000, 95);
fprintf(statsF,'95%%CI = [%.3f %.3f]; compare to 0.5, t(%i) = %.2f, p=%.4f\n', CI(1), CI(2), sts.df, sts.tstat, pval);

dists = allR.dualAccDistFromAllOrNone_PCs;
fprintf(statsF,'\nMean distance from serial model line = %.2f +/- %.2f\n\t', mean(dists), standardError(dists));
[~,pval,~,sts] = ttest(dists,0);
CI = boyntonBootstrap(@mean, dists, 5000, 95);
fprintf(statsF,'95%%CI = [%.3f %.3f]; compare to 0, t(%i) = %.2f, p=%.4f\n', CI(1), CI(2), sts.df, sts.tstat, pval);


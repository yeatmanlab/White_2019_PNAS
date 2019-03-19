%% function Fig1B_MeanBehavioralAOC() 
% Creates Figure 1B for White, Palmer, Boynton & Yeatman, PNAS 2019. 
% Plots the mean attention operating characteristic from proportion
% correct in semantic categorization task, collected in the scanner. 
% Requires the function plotAOCWithPredictions. 
% 
% by Alex L. White
% University of Washington, 2019

function Fig1B_MeanBehavioralAOC() 

%% set paths
%add whole analysis code directory to the path
analysisDir = fileparts(fileparts(which(mfilename)));
addpath(genpath(analysisDir));

p = getPaths(); 

%% load data 
resFile = fullfile(p.data,'AllSubjBehavior.mat'); 
load(resFile); 

%% extract data 
conds.cueCond = [3 2]; %focal, divided
conds.targSide = 2:3; %left, right
conds.congruent = 1; %not divided by congruency
conds.targPres = 1; %not divided by target presence/absence

as = squeeze(rAvg.PCs(conds.cueCond, conds.targSide, conds.congruent, conds.targPres)); 
es = squeeze(rAvg.SEM.PCs(conds.cueCond, conds.targSide, conds.congruent, conds.targPres));


%% plot 
figure; hold on;

fontSize = 10;
plotOpt.markSz = 7;

xColors = [46 63 153; 147 26 29]/255;
dualFillColor = [1 1 1];

plotOpt.edgeColors = xColors([1 1],:);
plotOpt.fillColors = [xColors(1,:); dualFillColor];

plotOpt.doAxLabels = true;
plotOpt.axLineWidth = 1;
plotOpt.datLineWidth = 1;

plotOpt.plotSerialPrediction = false;
plotOpt.axticks = 0.5:0.1:1; 
plotOpt.units = 'p(correct)';
plotOpt.doLegend = false;

for ati=1:length(plotOpt.axticks)
    if mod(ati,2)==1
        plotOpt.axtickLabels{ati} = '';
    else
        plotOpt.axtickLabels{ati} = sprintf('%.1f',plotOpt.axticks(ati));
    end
end

plotAOCWithPredictions(as,es,plotOpt);

set(gcf,'color','w','units','centimeters','pos',[5 5 5 5]);
exportfig(gcf,fullfile(p.figures,'Fig1B_AccuracyAOC.eps'),'Format','eps','bounds','loose','color','rgb','LockAxes',0,'FontMode','fixed','FontSize',fontSize);

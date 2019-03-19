%% Figure 4 for White, Palmer, Boynton & Yeatman, PNAS 2019
% Plots Attention Operating Characteristics (AOCs) from the across-subject
% mean BOLD responses in retinotopic cortex, and from channel responses in
% left VWFA-1 and VWFA-2. 
%
% by Alex L. White at the University of Washington, 2019

function Fig4_NeuronalAOCs()

close all;

%% set paths
%add whole analysis code directory to the path
analysisDir = fileparts(fileparts(which(mfilename)));
addpath(genpath(analysisDir));

p = getPaths(); 

%% open a file to print stats
statsFile = fullfile(p.stats,'Stats4_NeuronalAOCs.txt');
statsF = fopen(statsFile,'w');

fprintf(statsF,'Stats on the distributed-cue point in neuronal AOCs\n');
%% load mean BOLD data from all retinotopic areas

resFile = fullfile(p.results,'AllSubjMainExptResponses.mat');
load(resFile);

distribCueI = find(strcmp(allR.valsByIndex.condition,'Distributed cue'));
focalLeftI = find(strcmp(allR.valsByIndex.condition,'Focal cue left'));
focalRightI = find(strcmp(allR.valsByIndex.condition,'Focal cue right'));

areaNamesToPlot = {'V1','V2','V3','V4','VO','LO'};
[~,areasToPlot] = intersect(allR.valsByIndex.brainArea,areaNamesToPlot,'stable');
dats = allR.meanBetas;
dats = squeeze(dats(areasToPlot,:,:,:));

%AVERAGE OF RETINOTTOPIC AREAS
rs = squeeze(nanmean(dats,1));
areasAveraged = allR.valsByIndex.brainArea(areasToPlot);

nSubj = size(rs,ndims(rs));
%% plot options

fSize = 9;

xColor = [46 63 153]/255;
dualFillColor = [1 1 1];

plotOpt.markSz = 7;
plotOpt.edgeColors = xColor([1 1],:);
plotOpt.fillColors = [xColor(1,:); dualFillColor];
plotOpt.axLineWidth = 1;
plotOpt.datLineWidth = 1;
plotOpt.errorBarWidth = 1.5;
plotOpt.axlims = [0 0.12];
plotOpt.sideLabels = {'Right hemisphere','Left hemisphere'};

%% Make AOC from mean responses in retinotopic cortex, and report across-subject stats

%compute values for AOC, which for the focal-cue axes are
%the differences between attend-in and attend-out, and for the
%distributed-cue points are the differences between attend-both and
%attend-out
%the input to plotting function, "as", is 2x2 matrix:
%    rows = focal, distributed
%    columns = left side, right side
allAs = zeros(2,2,nSubj);
for hi=1:2
    if hi==1 %Left hemisphere
        attnInI = focalRightI;
        attnOutI = focalLeftI;
    else %rigth hemisphere
        attnInI = focalLeftI;
        attnOutI = focalRightI;
    end
    %focal cue points:
    allAs(1,hi,:) = rs(hi,attnInI,:)-rs(hi,attnOutI,:);
    %distributed cue points:
    allAs(2,hi,:) = rs(hi,distribCueI,:)-rs(hi,attnOutI,:);
end

%compute stats on distribute cue position for each subject
retinoDistsFromSwitchLine = NaN(1,nSubj);
retinoDistsFromUnlimited = NaN(1,nSubj);
goodSubj = false(1,nSubj);
for si=1:nSubj
    as = squeeze(allAs(:,:,si));
    %flip so that RH is on y-axis, LH on x-axis
    as = as(:,[2 1]);
    if ~any(isnan(as(:))) && ~any(as(1,:)<0)  %doesnt make sense to plot if some single-task values are negative!
        goodSubj(si) = true;
        [retinoDistsFromSwitchLine(si), ~, retinoDistsFromUnlimited(si)] = computeAOCPointDistances(as,[]);
    end
end

% stats
fprintf(statsF,'\n------------------------------------------------------\n');
fprintf(statsF,'ANALYSIS OF THE AVERAGE OF ALL RETINOTOPIC AREAS');
fprintf(statsF,'\n------------------------------------------------------\n');
fprintf(statsF,'Areas averaged:\n\t');
for avi=1:numel(areasAveraged)
    fprintf(statsF,'%s\t',areasAveraged{avi});
end
fprintf(statsF,'\n\n\n');

for distType = 1:2
    if distType==1
        dists = retinoDistsFromSwitchLine;
        fprintf(statsF,'(N that could be analyzed: %i)\n', sum(goodSubj));
        fprintf(statsF,'\nDistance from the nearest point on the serial switching line:\n');
        positiveDescription = '%i/%i above the line\n';
    elseif distType==2
        dists = retinoDistsFromUnlimited;
        fprintf(statsF,'\nDistance from the unlimited capacity point:\n');
        positiveDescription = '%i/%i better than that point, averaging over hemispheres\n';
    end
    
    distCI = boyntonBootstrap(@mean,dists(~isnan(dists)),5000,95);
    
    fprintf(statsF,'mean:\t %.3f\t SEM:\t %.3f\t 95%%CI=[%.3f %.3f]\n', ...
        nanmean(dists), standardError(dists), distCI(1), distCI(2));
    fprintf(statsF,positiveDescription, sum(dists>0), sum(goodSubj));
end

nGoodSubj = sum(all(all(~isnan(allAs),1),2));

as = nanmean(allAs,ndims(allAs));
es = nanstd(allAs,0,ndims(allAs))/sqrt(nGoodSubj);

%flip so that RH is on y-axis, LH on x-axis
as = as(:,[2 1]);
es = es(:,[2 1]);

%also compute distances in the averaged AOC
[meanRetinoDistFromSwitchLine, ~, meanRetinoDistFromUnlimited] = computeAOCPointDistances(as,[]);
fprintf(statsF,'\n\nIn the AVERAGE AOC from retinotopic cortex:\n\tDistance from the nearest point on the serial switching line = %.3f',meanRetinoDistFromSwitchLine);
fprintf(statsF,'\n\tDistance from the unlimited capacity point = %.3f\n',meanRetinoDistFromUnlimited);


%% Plot retinotopic AOC 
plotOpt.doAxLabels = true;
plotOpt.doLegend = false;
plotOpt.axlims = [0 0.12];
plotOpt.axticks = 0:0.03:0.12;
plotOpt.doXLabel = true;
plotOpt.doYLabel = true;

figure(1);
subplot(1,3,1); hold on;
plotAOCfromBOLDData(as,es,plotOpt);

title('Retinotopic areas');
        
set(gca,'LabelFontSizeMultiplier',1.0);
set(gca,'TitleFontSizeMultiplier',1.0);
set(gca,'TitleFontWeight','normal');

%% load channel respones for left VWFA-1 and VWFA-2

resFile = fullfile(p.results,'AllSubjChannelResponses.mat');
load(resFile);

areaNames = {'VWFA_1','VWFA_2'};
areaIs = zeros(size(areaNames));
for ai=1:length(areaNames)
    areaIs(ai) = find(strcmp(allR.valsByIndex.brainArea, areaNames{ai}));
    areaNames{ai}(areaNames{ai}=='_') = ' ';
end
nAreas = length(areaIs);

hemName = 'Left';
hemI = find(strcmp(allR.valsByIndex.hemisphere,hemName));
condIs = 1:3;
channIs = 1:2;

dats = squeeze(allR.twoChannelResponses(areaIs,hemI,channIs,condIs,:));

distribCueI = find(strcmp(allR.valsByIndex.condition,'Distributed cue'));
focalLeftI = find(strcmp(allR.valsByIndex.condition,'Focal cue left'));
focalRightI = find(strcmp(allR.valsByIndex.condition,'Focal cue right'));

%% compute VWFA channel response AOCs

plotOpt.axlims = [0 0.4];
plotOpt.axticks = 0:0.1:0.4;
plotOpt.sideLabels = {'Left channel','Right channel'};

for aii = 1:nAreas
    
    %pull out data for this brain area
    rs = squeeze(dats(aii, :, :, :));
    
    allAs = zeros(2,2,nSubj);
    for ki=1:2
        if ki==1 %Left channel
            attnInI = focalLeftI;
            attnOutI = focalRightI;
        else %right channel
            attnInI = focalRightI;
            attnOutI = focalLeftI;
        end
        %focal cue points:
        allAs(1,ki,:) = rs(ki,attnInI,:)-rs(ki,attnOutI,:);
        %distributed cue points:
        allAs(2,ki,:) = rs(ki,distribCueI,:)-rs(ki,attnOutI,:);
    end
    
    %compute stats on distribute cue position for each subject
    vwfaDistsFromSwitchLine = NaN(1,nSubj);
    vwfaDistsFromUnlimited = NaN(1,nSubj);
    goodSubj = false(1,nSubj);
    for si=1:nSubj
        as = squeeze(allAs(:,:,si));
        if ~any(isnan(as(:))) && ~any(as(1,:)<0) %doesnt make sense to plot if some single-task values are negative!
            goodSubj(si) = true;
            [vwfaDistsFromSwitchLine(si), ~, vwfaDistsFromUnlimited(si)] = computeAOCPointDistances(as,[]);
        end
    end
    
    % stats
    fprintf(statsF,'\n------------------------------------------------------\n');
    fprintf(statsF,'ANALYSIS OF THE AVERAGE OF %s %s', hemName, areaNames{aii});
    fprintf(statsF,'\n------------------------------------------------------\n');
    
    for distType = 1:2
        if distType==1
            dists = vwfaDistsFromSwitchLine;
            fprintf(statsF,'(N that could be analyzed: %i)\n', sum(goodSubj));
            fprintf(statsF,'\nDistance from the nearest point on the serial switching line:\n');
            positiveDescription = '%i/%i above the line\n';
        elseif distType==2
            dists = vwfaDistsFromUnlimited;
            fprintf(statsF,'\nDistance from unlimited capacity point:\n');
            positiveDescription = '%i/%i better than that point, averaging over hemispheres\n';
            
        end
        
        distCI = boyntonBootstrap(@mean,dists(~isnan(dists)),5000,95);
        
        fprintf(statsF,'mean:\t %.3f\t SEM:\t %.3f\t 95%%CI=[%.3f %.3f]\n', ...
            nanmean(dists), standardError(dists), distCI(1), distCI(2));
        
        fprintf(statsF,positiveDescription, sum(dists>0), sum(goodSubj));
    end
    
    nGoodSubj = sum(all(all(~isnan(allAs),1),2));
    
    as = nanmean(allAs,ndims(allAs));
    es = nanstd(allAs,0,ndims(allAs))/sqrt(nGoodSubj);
    
    
    %also compute distances in the averaged AOC
    [meanDistFromSwitchLine, ~, meanDistFromUnlimited] = computeAOCPointDistances(as,[]);
    
    fprintf(statsF,'\n\nIn the AVERAGE AOC:\n\tDistance from the nearest point on the serial switching line = %.3f',meanDistFromSwitchLine);
    fprintf(statsF,'\n\tDistance from the unlimited capacity point = %.3f\n',meanDistFromUnlimited);

    %% Plot channel response AOC
    figure(1);
    subplot(1,3,aii+1);
    
    plotOpt.axlims = [0 0.4];
    plotOpt.axticks = 0:0.1:0.4;
    plotOpt.doXLabel = true;
    plotOpt.doYLabel = true;
    
    plotAOCfromBOLDData(as,es,plotOpt);

    title([hemName ' ' areaNames{aii}]);
    set(gca,'LabelFontSizeMultiplier',1.0);
    set(gca,'TitleFontSizeMultiplier',1.0);
    set(gca,'TitleFontWeight','normal');
end

%save figure
set(gcf,'color','w','units','centimeters','pos',[5 5 12.4 4.6]);
figTtl = 'Fig4_MeanBOLDAOCs.eps';
exportfig(gcf,fullfile(p.figures,figTtl),'Format','eps','bounds','loose','color','rgb','LockAxes',0,'FontMode','fixed','FontSize',fSize);
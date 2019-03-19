%% function Stats2C_BOLDRespLME() 
% Printing statistics for White, Palmer, Boynton & Yeatman, PNAS 2019. 
% This function prints stats about mean BOLD responses in the main experiment, to assess 
% selective and divided spatial attention effects. 
% 
% by Alex L. White
% University of Washington, 2019

function Stats2C_BOLDRespLME() 
clc;

%% set paths
%add whole analysis code directory to the path
analysisDir = fileparts(fileparts(which(mfilename)));
addpath(genpath(analysisDir));

p = getPaths(); 

%% load data
resFile = fullfile(p.results,'AllSubjMainExptResponses.mat'); 
load(resFile); 


%% Open a file to print stats to:
statsFile = fullfile(p.stats,'Stats2C_MainExptBOLDStats.txt');
diary(statsFile);
statsF = fopen(statsFile,'w'); 

fprintf(1,'ANALYSIS OF BOLD response amplitudes, averaged over voxels in each ROI, as estimated by glmDenoise\n'); 
fprintf(1,'Unless otherwise stated, ''DummyVarCoding'' = ''effects'', so we can estimate main effects and interactions:\n');


%% extract data 
condIs = 1:3;
condLabs = allR.valsByIndex.condition(condIs);
nConds = numel(condLabs);

hemIs = 1:2;
hemLabs = allR.valsByIndex.hemisphere(hemIs);
nHems = length(hemIs);

%include all areas 
areasToPlot = 1:length(allR.valsByIndex.brainArea);
nAreas = length(areasToPlot);

%label the brain areas and take out underscores
areaLabs = allR.valsByIndex.brainArea(areasToPlot);
for ai=1:nAreas
    al = areaLabs{ai};
    al(al=='_') = '-';
    areaLabs{ai} = al;
end

nSubj = size(allR.meanBetas,ndims(allR.meanBetas));

bs = squeeze(allR.meanBetas(areasToPlot,hemIs,condIs,:));

%% print mean selective attention effects in each area 
fLI = find(strcmp(condLabs,'Focal cue left'));
fRI = find(strcmp(condLabs,'Focal cue right'));
dI = find(strcmp(condLabs,'Distributed cue'));

allSelectEs = NaN(nAreas, nHems, nSubj);
for hi=1:2
    if strcmp(hemLabs{hi},'Left')
        allSelectEs(:,hi,:) = bs(:,hi,fRI,:) - bs(:,hi,fLI,:); %focal cue contra - cue ipsi
    elseif strcmp(hemLabs{hi},'Right')
        allSelectEs(:,hi,:) = bs(:,hi,fLI,:) - bs(:,hi,fRI,:); %focal cue contra - cue ipsi
    end
end

%average over retinotopic areas
retinoNames = {'V1','V2','V3','V4','LO','VO'};
[~,areasToAvg] = intersect(areaLabs,retinoNames,'stable');
retinoSelectEs = mean(allSelectEs(areasToAvg,:,:),1);

[otherAreas,otherIs] = setdiff(areaLabs,retinoNames,'stable');
selectEs = cat(1,retinoSelectEs, allSelectEs(otherIs,:,:));

selectAreaNames = cat(1,{'Retinotopic'},otherAreas);

meanSelectEs = nanmean(selectEs,ndims(selectEs));
semSelectEs = standardError(selectEs,ndims(selectEs));

%also averaged over hems 
hemMeanSelectEs = squeeze(nanmean(selectEs,2)); 
meanHemMeanSelectEs = nanmean(hemMeanSelectEs, ndims(hemMeanSelectEs));
semHemMeanSelectEs = standardError(hemMeanSelectEs, ndims(hemMeanSelectEs));

fprintf(1,'\nMean (SEM) Selective Attention Effects (focal cue contralateral - ipsilateral) in each ROI:\n');
for ai=1:length(selectAreaNames)
    for hi=1:nHems
        fprintf(1,'\n%s %s: %.4f (%.4f)', hemLabs{hi}, selectAreaNames{ai}, meanSelectEs(ai,hi), semSelectEs(ai,hi));
    end
    if ai==1, fprintf(1,'\n\tAvged over hems: %.4f (%.4f)\n', meanHemMeanSelectEs(ai), semHemMeanSelectEs(ai)); end
end


%% print mean DIVIDED attention effects in each area 

allDivEs = NaN(nAreas, nHems, nSubj);
for hi=1:2
    if strcmp(hemLabs{hi},'Left')
        allDivEs(:,hi,:) = bs(:,hi,fRI,:) - bs(:,hi,dI,:); %focal cue contra - distributed
    elseif strcmp(hemLabs{hi},'Right')
        allDivEs(:,hi,:) = bs(:,hi,fLI,:) - bs(:,hi,dI,:); %focal cue contra - distributed
    end
end

%average over retinotopic areas
retinoNames = {'V1','V2','V3','V4','LO','VO'};
[~,areasToAvg] = intersect(areaLabs,retinoNames,'stable');
retinoDivEs = mean(allDivEs(areasToAvg,:,:),1);

[otherAreas,otherIs] = setdiff(areaLabs,retinoNames,'stable');
dividEs = cat(1,retinoDivEs, allDivEs(otherIs,:,:));

selectAreaNames = cat(1,{'Retinotopic'},otherAreas);

meanDividEs = nanmean(dividEs,ndims(dividEs));
semDividEs = standardError(dividEs,ndims(dividEs));

%also averaged over hems 
hemMeanDividEs = squeeze(nanmean(dividEs,2)); 
meanHemMeanDividEs = nanmean(hemMeanDividEs, ndims(hemMeanDividEs));
semHemMeanDividEs = standardError(hemMeanDividEs, ndims(hemMeanDividEs));

fprintf(1,'\n\nMean (SEM) DIVIDED Attention Effects (focal cue contralateral - distributed cue) in each ROI:\n');
for ai=1:length(selectAreaNames)
    for hi=1:nHems
        fprintf(1,'\n%s %s: %.4f (%.4f)', hemLabs{hi}, selectAreaNames{ai}, meanDividEs(ai,hi), semDividEs(ai,hi));
    end
    if ai==1,    fprintf(1,'\n\tAvged over hems: %.4f (%.4f)\n', meanHemMeanDividEs(ai), semHemMeanDividEs(ai)); end
end

%% create a table
nRows       = numel(bs);
subject    = NaN(nRows,1); 
region     = cell(nRows,1); 
hemisphere = cell(nRows,1); 
cue        = cell(nRows,1); 
bold       = NaN(nRows,1);

r = 0;
for si = 1:nSubj
    for ai = 1:nAreas
        for hi = 1:nHems
            for ci = 1:nConds
                r = r+1;
                subject(r) = si;

                region(r) = areaLabs(ai);
                hemisphere(r) = hemLabs(hi);
                
                if strcmp(condLabs{ci},'Distributed cue')
                    cue(r) = {'distributedCue'};
                else
                    if strcmp(hemisphere{r},'Left')
                        if strcmp(condLabs{ci}, 'Focal cue left')
                            cue(r) = {'focalIpsilateral'};
                        elseif strcmp(condLabs{ci}, 'Focal cue right')
                            cue(r) = {'focalContralateral'};
                        end
                    elseif strcmp(hemisphere{r},'Right')
                        if strcmp(condLabs{ci}, 'Focal cue right')
                            cue(r) = {'focalIpsilateral'};
                        elseif strcmp(condLabs{ci}, 'Focal cue left')
                            cue(r) = {'focalContralateral'};
                        end
                    end
                end
                
                bold(r) = bs(ai,hi,ci,si);
                
            end
        end
    end
end

subject = categorical(subject);
d = table(subject, region, hemisphere, cue, bold);


fprintf(1,'\n\n------------------------------------------------------\n');
fprintf(1,'\n------------------------------------------------------\n');
fprintf(1,'ANALYSIS OF BETA WEIGHTS FROM INDIVIDUAL HEMISPHERES\n');
fprintf(1,'Coding the cue as distributed, focal contralateral, or focal ipsilateral\n');
fprintf(1,'Averaging over retinotopic areas');
fprintf(1,'\n------------------------------------------------------\n');

%% Fit LME
%1. Three-way effects of cue, hemisphere, and region
fprintf(1,'\n------------------------------------------------------\n');
fprintf(1,'Linear Mixed Effects model fit for three-way effects of region (3), hemisphere (2), and cue (3):\n');

eqtn = 'bold ~ region * hemisphere * cue + (1 | subject)';
subst = true(size(d,1),1); %all data

lme1 = fitlme(d(subst,:),eqtn,'DummyVarCoding','effects');
lme1.disp;
fprintf(1,'\n---Corresponding ANOVA---\n'); 
display(lme1.anova);

%2. Three-way effects of cue, hemisphere, and region, excluding dual-task
fprintf(1,'\n------------------------------------------------------\n');
fprintf(1,'Linear Mixed Effects model fit for three-way effects of region (8), hemisphere (2), and cue (2), EXCLUDING distributed-cue condition:\n');

eqtn = 'bold ~ region * hemisphere * cue + (1 | subject)';
subst = ~strcmp(d.cue,'distributedCue');

lme2 = fitlme(d(subst,:),eqtn,'DummyVarCoding','effects');
display(lme2);

fprintf(1,'\n---Corresponding ANOVA---\n'); 
display(lme2.anova);

%do it again but with random cue effects by subject 
fprintf(1,'\n------------------------------------------------------\n');
fprintf(1,'Same model, but with cue effect also varying by subject:\n');

eqtn3 = 'bold ~ region * hemisphere * cue + (cue | subject)';

lme3 = fitlme(d(subst,:),eqtn3,'DummyVarCoding','effects');
display(lme3);

fprintf(1,'\n---Corresponding ANOVA---\n'); 
display(lme3.anova);
fprintf(1,'Model comparison to test if that is necessary:\n');
compare(lme2,lme3)

%2b. Three-way effects of cue, hemisphere, and region, just focal
%contralateral vs distributed (excluding focal ipsi)
fprintf(1,'\n------------------------------------------------------\n');
fprintf(1,'Linear Mixed Effects model fit for three-way effects of region (3), hemisphere (2), and cue (2), just focal contralateral vs distributed (excluding focal ipsi):\n');

eqtn = 'bold ~ region * hemisphere * cue + (1 | subject)';
subst = ~strcmp(d.cue,'focalIpsilateral');

lme2b = fitlme(d(subst,:),eqtn,'DummyVarCoding','effects');
display(lme2b);

fprintf(1,'\n---Corresponding ANOVA---\n'); 
display(lme2b.anova);


%4. Three-way effects of cue, hemisphere, and region, only in VWFA-1 and
%VWFA_2
fprintf(1,'\n------------------------------------------------------\n');
fprintf(1,'Linear Mixed Effects model fit for three-way effects of region (2), hemisphere (2), and cue (3), only in VWFA-1 and VWFA_2:\n');

eqtn = 'bold ~ region * hemisphere * cue + (1 | subject)';
subst = strcmp(d.region,'VWFA-1') | strcmp(d.region,'VWFA-2');

lme4 = fitlme(d(subst,:),eqtn,'DummyVarCoding','effects');
display(lme4);

fprintf(1,'\n---Corresponding ANOVA---\n'); 
display(lme4.anova);

%5. Three-way effects of cue, hemisphere, and region, only in VWFA-1 and
%VWFA_2, excluding dual-task 
fprintf(1,'\n------------------------------------------------------\n');
fprintf(1,'Linear Mixed Effects model fit for three-way effects of region (2), hemisphere (2), and cue (2), only in VWFA-1 and VWFA_2, excluding distributed cue\n');

eqtn = 'bold ~ region * hemisphere * cue + (1 | subject)';
subst = (strcmp(d.region,'VWFA-1') | strcmp(d.region,'VWFA-2')) & ~strcmp(d.cue,'distributedCue');

lme5 = fitlme(d(subst,:),eqtn,'DummyVarCoding','effects');
display(lme5);

fprintf(1,'\n---Corresponding ANOVA---\n'); 
display(lme5.anova);


%6. Three-way effects of cue, hemisphere, and region, only in VWFA-1 and
%VWFA_2, excluding single-task ipsi (divided attn effects)
fprintf(1,'\n------------------------------------------------------\n');
fprintf(1,'Linear Mixed Effects model fit for three-way effects of region (2), hemisphere (2), and cue (2), only in VWFA-1 and VWFA_2, just focal contralateral vs distributed (excluding focal ipsi)\n');

eqtn = 'bold ~ region * hemisphere * cue + (1 | subject)';
subst = (strcmp(d.region,'VWFA-1') | strcmp(d.region,'VWFA-2')) & ~strcmp(d.cue,'focalIpsilateral');

lme6 = fitlme(d(subst,:),eqtn,'DummyVarCoding','effects');
display(lme6);

fprintf(1,'\n---Corresponding ANOVA---\n'); 
display(lme6.anova);

diary off;

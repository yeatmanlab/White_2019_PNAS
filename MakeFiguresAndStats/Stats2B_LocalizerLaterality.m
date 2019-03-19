%% function Stats2B_LocalizerLaterality() 
% Printing statistics for White, Palmer, Boynton & Yeatman, PNAS 2019. 
% This function prints stats about mean BOLD responses to single words on
% the left and right of fixation from the localizer scans, in order to
% assess the hemifield selectivity of each brain region. 
% 
% by Alex L. White
% University of Washington, 2019

function Stats2B_LocalizerLaterality()

%% set paths
%add whole analysis code directory to the path
analysisDir = fileparts(fileparts(which(mfilename)));
addpath(genpath(analysisDir));

p = getPaths(); 

%% load data

resFile = fullfile(p.results,'AllSubjLocalizerResponses.mat');
load(resFile);

%% open stats folder 
statsFileName = fullfile(p.stats,'Stats2B_VoxelLaterality.txt');
statsF = fopen(statsFileName,'w');
diary(statsFileName);

fprintf(statsF,'STATISTICS ON VOXEL HEMIFIELD SELECTIVITY\n\n');
fprintf(statsF,'Comparing responses to words on the left and right sides of fixation, in VWFA_1 and VWFA_2\n');
fprintf(statsF,'\nUsing the lateralization index: LI = 1 - ipsilateral/contralateral\n');

%% extract data
areasToPlot = [find(strcmp(allR.valsByIndex.brainArea,'VWFA_1')) find(strcmp(allR.valsByIndex.brainArea,'VWFA_2'))];
nAreas = length(areasToPlot);

hemIs = 1:2;
hemLabs = allR.valsByIndex.hemisphere(hemIs);
nHems = length(hemIs);

condsToPlot = {'left_word','right_word'};

nConds = numel(condsToPlot);
condIs = zeros(1,nConds);
condLabs = cell(1,nConds);
for ci=1:nConds
    condIs(ci) = find(strcmp(allR.valsByIndex.condition,condsToPlot{ci}));
    rl = condsToPlot{ci};
    rl(rl=='_') = ' ';
    condLabs{ci} = rl;
end

%label the brain areas and take out hyphens
areaLabs = allR.valsByIndex.brainArea(areasToPlot);
for ai=1:nAreas
    al = areaLabs{ai};
    al(al=='_') = '-';
    areaLabs{ai} = al;
end

nSubj = size(allR.meanBetas,ndims(allR.meanBetas));

%% compute lateralization index: 1-(ipsilateral/contralateral)
allbs = squeeze(allR.meanBetas(areasToPlot,hemIs,condIs,:));
LIs = NaN(nAreas, nHems, nSubj);
for hi=1:2
    % for LH, the ratio is right/left
    % for RH, the ratio is left/right
    LIs(:,hi,:) = 1 - (allbs(:,hi,hi,:)./allbs(:,hi,3-hi,:));
end

meanLIs = nanmean(LIs,ndims(LIs));
semLIs = standardError(LIs,ndims(LIs));

%print mean LIs in each region 
fprintf(1,'\nMean (SEM) Lateralization Indices:');
for hii=1:nHems
    for aii=1:nAreas
        fprintf(1,'\n%s %s: %.2f (%.2f)', hemLabs{hii}, areaLabs{aii}, meanLIs(aii,hii), semLIs(aii,hii));
    end
end

%% put those into a table
nRows = nAreas*nHems*nSubj;

LI = NaN(nRows,1);
subject = NaN(nRows,1);
region = cell(nRows,1);
hemisphere = cell(nRows,1);

rowI = 0;
for aii=1:nAreas
    for hii = 1:nHems
        for sii=1:nSubj
            rowI = rowI+1;
            LI(rowI) = LIs(aii, hii, sii);
            subject(rowI) = sii;
            region{rowI} = areaLabs{aii};
            hemisphere{rowI} = hemLabs{hii};
        end
    end
end

subject = categorical(subject);
T = table(LI, region, hemisphere, subject);

%% Next: run LME

fprintf(1,'\n\n\nLME for effects of region*hemisphere on LI\n');
eqtn = 'LI ~ region * hemisphere + (1 | subject)';

lme = fitlme(T,eqtn,'DummyVarCoding','effects');
lme.disp;
fprintf(1,'\n---Corresponding ANOVA---\n');
display(lme.anova);

fprintf(1,'\n\n\nLME for effects of region on LI in just LEFT hemisphere VWFA_1 vs VWFA_2\n');

eqtn = 'LI ~ region + (1 | subject)';
subst = strcmp(T.hemisphere,'Left');
lme = fitlme(T(subst,:),eqtn,'DummyVarCoding','effects');
lme.disp;
fprintf(1,'\n---Corresponding ANOVA---\n');
display(lme.anova);

diary off;

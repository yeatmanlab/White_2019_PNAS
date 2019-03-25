%% function ComputeMeanBOLDResps() 
% 
% Analyzes data for White, Palmer, Boynton & Yeatman, PNAS 2019
% This function loads in a table (T) that contains the response of each
% individual voxel in each condition of the localizer scans and the main
% experiment. 
% This function then computes the across-voxel mean responses in each ROI, in each subject, in the
% localizer scans and main experimental scans. It also computes mean
% cross-validated r2  and snr statistics from glmDenoise. 
% For the localizer scans, the results are accumulated into a single variable "allR" that is saved
% into a mat file "AllSubjLocalizerResponses.mat" in the "data" folder. 
% 
% For the main scans, the results are accumulated into a single variable "allR" that is saved
% into a mat file "AllSubjMainExptResponses.mat" in the "data" folder. 
% 

% by Alex L. White at the University of Washington, 2019

function ComputeMeanBOLDResps() 

%% set paths
%add whole analysis code directory to the path
analysisDir = fileparts(fileparts(which(mfilename)));
addpath(genpath(analysisDir));

p = getPaths(); 


%% load big table of individual voxel responses
tableFileName = fullfile(p.data,'AllSubjectVoxelResponseTable.csv'); 
T = readtable(tableFileName);

brainAreas = unique(T.region); 
hemispheres = unique(T.hemisphere); 

%% analyze localizer data 

%which conditions to extract 
stimConds = {'resp_wordL','resp_wordR'}; %also includes resp_scramL and resp_scramR 

%alternate labels
stimCondLabels = cell(size(stimConds)); 
stimCondLabels(strcmp(stimConds,'resp_wordL')) = {'left_word'};
stimCondLabels(strcmp(stimConds,'resp_wordR')) = {'right_word'};

nAreas = length(brainAreas); 
nHems = length(hemispheres); 
nConds = length(stimConds);
nSubj = max(T.subject); 

meanBetas = NaN(nAreas, nHems, nConds, nSubj); 
numVoxels = NaN(nAreas, nHems, nSubj); 
xvalRSqr  = NaN(nAreas, nHems, nSubj); 
snr       = NaN(nAreas, nHems, nSubj); 

for ai=1:nAreas
    for hi = 1:nHems
        for si = 1:nSubj
            voxelRows = strcmp(T.region,brainAreas{ai}) & strcmp(T.hemisphere,hemispheres{hi}) & T.subject==si;
            
            %count voxels in this ROI for this subject
            numVoxels(ai, hi, si) = sum(voxelRows);
            
            if sum(voxelRows)>0
                
                %cross-validated R2 from glmDenoise, averaged across voxels
                xvalRSqr(ai, hi, si)  = mean(T.xvalRSqr_locz(voxelRows));
                
                %SNR from glmDenosie, averaged across voxels
                snr(ai, hi, si)       =  mean(T.snrAfter_locz(voxelRows));
                
                %beta weights in each stimulus condiiton, averaged across voxels
                for ci = 1:nConds
                    %extract responses for this condition
                    eval(sprintf('betas = T.%s(voxelRows);', stimConds{ci}));
                    %average over voxels
                    meanBetas(ai, hi, ci, si) = mean(betas);
                end
            end
        end
    end
end
           
%create the 'valsByIndex' structure, which labels each level of each
%dimension of the meanBetas matrix. 
allR.valsByIndex.brainArea = brainAreas; 
allR.valsByIndex.hemisphere = hemispheres; 
allR.valsByIndex.condition = stimCondLabels;

%put the rest of the data in allR 
allR.meanBetas = meanBetas;
allR.numVoxels = numVoxels; 
allR.xvalRSqr = xvalRSqr; 
allR.snr = snr;

%save this file 
resFileName = fullfile(p.results,'AllSubjLocalizerResponses.mat'); 
save(resFileName,'allR');

%% analyze main experiment data 

clear allR; 

%conditions to extract from the table 
cueConds = {'resp_distributedCue','resp_focalCueLeft','resp_focalCueRight'};

%alternate labels
cueCondLabels = cell(size(cueConds)); 
cueCondLabels(strcmp(cueConds,'resp_distributedCue')) = {'Distributed cue'};
cueCondLabels(strcmp(cueConds,'resp_focalCueLeft')) = {'Focal cue left'};
cueCondLabels(strcmp(cueConds,'resp_focalCueRight')) = {'Focal cue right'};

nAreas = length(brainAreas); 
nHems = length(hemispheres); 
nConds = length(cueConds);
nSubj = max(T.subject); 

meanBetas = NaN(nAreas, nHems, nConds, nSubj); 
numVoxels = NaN(nAreas, nHems, nSubj); 
xvalRSqr  = NaN(nAreas, nHems, nSubj); 
snr       = NaN(nAreas, nHems, nSubj); 

for ai=1:nAreas
    for hi = 1:nHems
        for si = 1:nSubj
            voxelRows = strcmp(T.region,brainAreas{ai}) & strcmp(T.hemisphere,hemispheres{hi}) & T.subject==si;
            
            %count voxels in this ROI for this subject
            numVoxels(ai, hi, si) = sum(voxelRows);
            
            if sum(voxelRows)>0
                
                %cross-validated R2 from glmDenoise, averaged across voxels
                xvalRSqr(ai, hi, si)  = mean(T.xvalRSqr_main(voxelRows));
                
                %SNR from glmDenosie, averaged across voxels
                snr(ai, hi, si)       =  mean(T.snrAfter_main(voxelRows));
                
                %beta weights in each stimulus condiiton, averaged across voxels
                for ci = 1:nConds
                    %extract responses for this condition
                    eval(sprintf('betas = T.%s(voxelRows);', cueConds{ci}));
                    %average over voxels
                    meanBetas(ai, hi, ci, si) = mean(betas);
                end
            end
        end
    end
end

%create the 'valsByIndex' structure, which labels each level of each
%dimension of the meanBetas matrix. 
allR.valsByIndex.brainArea = brainAreas; 
allR.valsByIndex.hemisphere = hemispheres; 
allR.valsByIndex.condition = cueCondLabels;

allR.meanBetas = meanBetas;
allR.numVoxels = numVoxels; 
allR.xvalRSqr = xvalRSqr; 
allR.snr = snr;

%save this file 
resFileName = fullfile(p.results,'AllSubjMainExptResponses.mat'); 
save(resFileName,'allR');

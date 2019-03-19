%% function ComputeMeanBOLDResps_WordFreq()
%
% Analyzes data for White, Palmer, Boynton & Yeatman, PNAS 2019
% This function loads in a table (T) that contains the response of each
% individual voxel in each condition of the main experiment, from a GLM
% that divided trials according to cue, the lexical frequency bin of the
% left word (low or high) and the lexical frequency bin of the right word
% (low or high). 
% This function then computes the across-voxel mean responses in each ROI, in each subject
% It also computes mean cross-validated r2  and snr statistics from glmDenoise. 
% The results are accumulated into a single variable "allR" that is saved
% into a mat file "AllSubjWordFreqResponses.mat" in the "data" folder. 
%
% by Alex L. White, University of Washington, 2019

function ComputeMeanBOLDResps_WordFreq()

%% set paths
%add whole analysis code directory to the path
analysisDir = fileparts(fileparts(which(mfilename)));
addpath(genpath(analysisDir));

p = getPaths(); 
%% load big table 
tableFileName = fullfile(p.data,'AllSubjectVoxelResponseTable_WordFreq.mat'); 
load(tableFileName,'T');

brainAreas = unique(T.region); 
hemispheres = unique(T.hemisphere); 


%% analyze main experiment data 

clear allR; 

%condition names to extract from the table 
conditionNames = {'distrib_leftLow_rightLow';
    'distrib_leftLow_rightHigh';
    'distrib_leftHigh_rightLow';
    'distrib_leftHigh_rightHigh';
    'focalLeft_leftLow_rightLow';
    'focalLeft_leftLow_rightHigh';
    'focalLeft_leftHigh_rightLow';
    'focalLeft_leftHigh_rightHigh';
    'focalRight_leftLow_rightLow';
    'focalRight_leftLow_rightHigh';
    'focalRight_leftHigh_rightLow';
    'focalRight_leftHigh_rightHigh'};


nAreas = length(brainAreas); 
nHems = length(hemispheres); 
nConds = length(conditionNames);
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
                    eval(sprintf('betas = T.resp_%s(voxelRows);', conditionNames{ci}));
                    %average over voxels
                    meanBetas(ai, hi, ci, si) = mean(betas);
                end
            end
        end
    end
end
     
allR.valsByIndex.brainArea = brainAreas; 
allR.valsByIndex.hemisphere = hemispheres; 
allR.valsByIndex.condition = conditionNames;

allR.meanBetas = meanBetas;
allR.numVoxels = numVoxels; 
allR.xvalRSqr = xvalRSqr; 
allR.snr = snr;

%save this file 
resFileName = fullfile(p.results,'AllSubjWordFreqResponses.mat'); 
save(resFileName,'allR');

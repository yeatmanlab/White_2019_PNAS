%% function SpatialEncodingModel_WordFreq() 
% Analyzes data for White, Palmer, Boynton & Yeatman, PNAS 2019
%
% This function run the spatial encoding model on each ROI in each subject,
% dividing trials by cue condition, the lexical frequency bin of the
% left word (low or high) and the lexical frequency bin of the right word
% (low or high). 
%
% First it loads in a big table, T, which has responses of each voxel in
% each condition of the localizer and each condition of the main attention
% experiment. 
%
% Then for each ROI, in each subject, it fits a linear model that assumes
% the responses in the main experiment are a weighted sum of responses in
% the localizer. 
%
% The two-channel model assumes that each voxels response is a weighted sum
% of two spatial channels, each of which was maximally and selecitvely
% stimulated by words at one location in the localizer scan. 
%
% The one-channel model assumes that each voxel's response is just a scaled
% version of the mean of its responses to left and right words in the
% localizer scan. 
%
% In addition to estimating channel responses for those two models in each
% cue condition, it also computes R2 for each model, and an adjusted R2
% that takes into account the number of free parameters. 
%
% All results are collected in a structure allR that is saved as
% AllSubjChannelResponses_WordFreq.mat, in the 'data' folder.
% 
% by Alex L. White, University of Washington, 2019
%
function SpatialEncodingModel_WordFreq() 

%% set paths
%add whole analysis code directory to the path
analysisDir = fileparts(fileparts(which(mfilename)));
addpath(genpath(analysisDir));

p = getPaths(); 

%% load big table 
tableFileName = fullfile(p.data,'AllSubjectVoxelResponseTable_WordFreq.mat'); 
load(tableFileName,'T');

%% select data
brainAreas = {'VWFA_1','VWFA_2'}; 
hemispheres = unique(T.hemisphere); 

nAreas    = length(brainAreas); 
nHems     = length(hemispheres); 
nChannels = 2;
nSubj = max(T.subject); 

loczCondsToExtract = {'resp_wordL','resp_wordR'};
channelLabs = cell(1,nChannels); 
channelLabs(strcmp(loczCondsToExtract, 'resp_wordL')) = {'Left'}; 
channelLabs(strcmp(loczCondsToExtract, 'resp_wordR')) = {'Right'}; 

mainCondsToExtract = {'distrib_leftLow_rightLow';
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

nAttnConds = length(mainCondsToExtract); 

%% run inverted encoding model for all subjects and ROIs
twoChannelResponses = NaN(nAreas, nHems, nChannels, nAttnConds, nSubj);
oneChannelResponses = NaN(nAreas,nHems, nAttnConds, nSubj);

nModelTypes = 2; %1-channel, 2-channel
nRSqrTypes = 2; %1=regular, 2=adjusted 

rSqrs = NaN(nAreas, nHems, nModelTypes, nRSqrTypes, nSubj); 

for si=1:nSubj
    for ai=1:nAreas
        for hi=1:nHems
            voxelRows = strcmp(T.region,brainAreas{ai}) & strcmp(T.hemisphere,hemispheres{hi}) & T.subject==si;
             
            if sum(voxelRows)>0 %only run model if this subject has this ROI 
                
                %voxel weights: responses to left and right words in the localizer scan
                W = NaN(sum(voxelRows), nChannels);
                for ki=1:nChannels
                    eval(sprintf('W(:,ki) = T.%s(voxelRows);', loczCondsToExtract{ki}));
                end
                
                %Voxel responses in each cue condition of the main experiment
                D = NaN(sum(voxelRows), nAttnConds);
                for ci=1:nAttnConds
                    eval(sprintf('D(:,ci) = T.resp_%s(voxelRows);', mainCondsToExtract{ci}));
                end
                
                %fit the one-channel model
                modelType = 1;
                [C, rSqr, rSqrAdj] = linearRegressionWithStats(mean(W,2), D);
                oneChannelResponses(ai,hi,:,si) = C;
                rSqrs(ai,hi,modelType,1,si) = rSqr;
                rSqrs(ai,hi,modelType,2,si) = rSqrAdj;
                
                %fit the two-channel model
                modelType = 2;
                [C, rSqr, rSqrAdj] = linearRegressionWithStats(W, D);
                twoChannelResponses(ai,hi,:,:,si) = C;
                rSqrs(ai,hi,modelType,1,si) = rSqr;
                rSqrs(ai,hi,modelType,2,si) = rSqrAdj;
            end
        end
    end
end

allR.oneChannelResponses = oneChannelResponses;
allR.twoChannelResponses = twoChannelResponses;

allR.valsByIndex.brainArea  = brainAreas; 
allR.valsByIndex.hemisphere = hemispheres; 
allR.valsByIndex.channel    = channelLabs;

allR.valsByIndex.condition = mainCondsToExtract;

allR.rSqrs = rSqrs;
allR.rSqrValsByIndex.brainArea = brainAreas;
allR.rSqrValsByIndex.hemisphere = hemispheres; 
allR.rSqrValsByIndex.modelNChannels = [1 2]; 
allR.rSqrValsByIndex.rSqrType  = {'regular','adjusted'};

%% save
%save this file 
resFileName = fullfile(p.results,'AllSubjChannelResponses_WordFreq.mat'); 
save(resFileName,'allR');



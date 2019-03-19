%% function StatsOnVoxelWeightCorrelations
% Printing statistics for White, Palmer, Boynton & Yeatman, PNAS 2019. 
% This function analyzes the correlations between single voxel responses to
% words on the left and right. It also examines the standard deviation of the differences 
% between respones to words on the left and right. 
% This is a control analysis to address the concern that left VWFA-2
% appears to lack 2 channels only because of lower SNR of our measurements.
% 
% % by Alex L. White
% University of Washington, 2019

% In more detail: 
% For each ROI in each subject, take the correlation across voxels of
% response to single word (L) on the left vs response to single words on the
% right (R). 
% 
% The one-channel model is that L = aR + noise
% The two-channel model is that L = aR + reliableDifference + noise
%    where "a" is some scalar difference reflecting mean hemifield asymmetry 
%
% So corr(L,R) must be lower for the two-channel model than the one-channel
% model
% 
function StatsOnVoxelWeightCorrelations()

%% set paths
%add whole analysis code directory to the path
analysisDir = fileparts(fileparts(which(mfilename)));
addpath(genpath(analysisDir));

p = getPaths(); 

%% load big table
tableFileName = fullfile(p.data,'AllSubjectVoxelResponseTable.mat');
load(tableFileName,'T');

%% stats file
statsF = fopen(fullfile(p.stats,'VoxelWeightCorrelations.txt'),'w');

fprintf(statsF,'ACROSS-VOXEL CORRELATIONS BETWEEN RESPONSES TO SINGLE WORDS ON THE LEFT AND SINGLE WORDSON THE RIGHT\n');
fprintf(statsF,'Computed on each subject separately, then averaged across subjects');
fprintf(statsF,'\nAlso examining the across-voxel standard deviation of *differences* in response to words on the left and words on the right\n');

%% extract voxels and compute stats

areaNames = {'VWFA_1','VWFA_2'};
hemNames = {'Left','Right'};

nAreas = length(areaNames);
nHems = length(hemNames);

uSubjs = unique(T.subject);
nSubj = numel(uSubjs);
for ai=1:nAreas
    for hi=1:nHems
        
        corrRhos = NaN(1,nSubj); 
        meanLRDiffs = NaN(1,nSubj);
        stdLRDiffs = NaN(1,nSubj);         
        
        for si=1:nSubj
            subst = strcmp(T.region,areaNames{ai}) & strcmp(T.hemisphere, hemNames{hi}) & T.subject==uSubjs(si);
            if sum(subst)>0
                
                %first compute correlation between responses to words on left and words on right
                wL = T.resp_wordL(subst,:);
                wR = T.resp_wordR(subst,:);
                
                corrRhos(si) = corr(wR,wL);
                
                %then stats on the differences between responses to words on the left and words on the right 
                lrDiffs = wL - wR; 
                meanLRDiffs(si) = mean(lrDiffs); 
                stdLRDiffs(si)  = std(lrDiffs); 
            end
        end
        
        fprintf(statsF,'\nFor %s %s:\n',areaNames{ai},hemNames{hi});
        fprintf(statsF,'\tMean (sem) rho = %.4f (%.3f)\n', nanmean(corrRhos), standardError(corrRhos));
        
        fprintf(statsF,'\n\tMean of across-voxel mean differences between wL and wR: %.4f (%.3f)', nanmean(meanLRDiffs), standardError(meanLRDiffs));
        fprintf(statsF,'\n\tMean of across-voxel standard deviation of differences between wL and wR: %.4f (%.3f)\n', nanmean(stdLRDiffs), standardError(stdLRDiffs));
    end
end



%% Script to process raw data from 
%% White, Palmer, Boynton & Yeatman, 2019 
% 
% Alex L. White, University of Washington 2019 
% 
% The functions used below load in tables of pre-processed individual voxel responses. 
% Those are stored in the "data" folder. The functions also assemble the
% results into mat files that are stored in the "results" folder. 
% Note that the behavioral data come already analysed, in mat files stored
% in the same "data" folder: AllSubjBehavior.mat and AllSubjBehavior_WordFreq.mat

%% 1. Compute mean BOLD responses in each ROI 
ComputeMeanBOLDResps; 

%% 2. Run the spatial encoding model 
SpatialEncodingModel; 

%% 3. For the analysis of word frequency effects, compute mean BOLD responses
ComputeMeanBOLDResps_WordFreq; 

%% 4. For the analysis of word frequency effects, run the spatial encoding model 
SpatialEncodingModel_WordFreq; 


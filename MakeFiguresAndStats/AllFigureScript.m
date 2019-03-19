%% Script to make all figures and print all statistics
%% for White, Palmer, Boynton & Yeatman, PNAS 2019 
% 
% by Alex L. White, University of Washington, 2019 
% 
% This script calls functions that load in results files from the "results"
% folder, save figures to the "Figures" folder, and print statistics to the
% "stats" folder. 
%%

Fig1B_MeanBehavioralAOC;
Stats1B_BehavioralAccuracy;

Fig2B_HemifieldSelectivityBarPlot;
Stats2B_LocalizerLaterality;

Fig2C_BOLDBarPlot;
Stats2C_BOLDRespLME;

Fig3_LeftHemChannelResps; 
Stats3_ChannelRespLME;

StatsOnVoxelWeightCorrelations;

Fig4_NeuronalAOCs;

Fig5_LexicalFrequency;

FigS1_LocalizerResps_2DHist;

FigS2_BOLDBarPlot_Retinotopic;

FigS3_EncodingModelAdjustedR2s;

FigS4_RightHemChannelResps;

FigS5_VoxelSpatialVsAttnSelectivity;

FigS6_LVWFA1_ChannelFreqEffects;

printROINVoxelTable;
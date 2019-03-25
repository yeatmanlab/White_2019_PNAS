# [Parallel spatial channels converge at a bottleneck in anterior word-selective cortex](https://www.biorxiv.org/content/10.1101/508846v1)
## White AL, Palmer J, Boynton GM, Yeatman JD

Code and data to reproduce figures and statistics from White, Palmer, Boynton and Yeatman (PNAS 2019).

Brief instructions: 
To produce the figures, you must run two Matlab scripts: 
1. AllAnalysisScript, which processes the raw data and saves aggregate results. 
2. AllFigureScript, which makes all the figures and also prints statistics to txt files. 

Detailed information on the contents of this repository: 

This repo contains 4 subfolders: 

(1) AnalyzeData: Matlab code to process the fMRI data. AllAnalysisScript.m calls four functions (also within that folder) to analyze different aspects of the data. Each function loads in a table T from the Data folder, which has each individual voxel's response in each condition of the localizer and main attention experiments. The results are saved in .mat files in a Results folder, which is created in the same parent folder as AnalyzeData/. 

Functions called by AllAnalysisScript: 
    ComputeMeanBOLDResps.m 
    SpatialEncodingModel.m
    ComputeMeanBOLDResps_WordFreq.m
    SpatialEncodingModel_WordFreq.m

More information is in the header of each .m file. 

(2) Data: .mat and .csv files containing data from each participant. There are eight files in this folder: 

PC.csv: a table of proportion correct for each participant in each main condition of the experiment, collected in the scanner. There are four columns: DistributedCue_Left, DistributedCue_Right, FocalCue_Left, FocalCue_Right. There is one row for each participant, and they are ordered in the same way as in the MRI data tables. 

dprimes.csv: a table of d' for each subject, formatted like PCs.csv. To compute d', we treated "living" words as the "signal", such that hits are when the participant correctly reports "living" and false alarms are when the participant incorrectly reports "living". Then d' is computed in the standard way from hit and false alarm rates in each condiiton (d' = norminv(HR) - norminv(FAR)). 

betas.csv: a table of beta values formatted like PCs.csv. In the signal detection framework, beta is the participant's criterion for reporting "target present" expressed as the likelihood ratio p(target present)/p(target absent). Again, "targets" are words from the "living" category. A beta of 1 is neutral, optimal given the assumption of equal-variance Gaussian distributions of evidence. Betas greater than 1 are conservative. 

meanCorrRTs: a table of geometric mean reaction times on correct trails, formatted like PCs.csv. The geometric mean is computed by first taking the log (base 10) of the reaction time on all correct trials, then averaging, then inverting the logarithm. The reaction time is in seconds, relative to the onset of the post-cue.

AllSubjBehavior.mat: 
    This mat file contains more information about task performance collected in the scanner. This file contains two variables: allR and rAvg. The same information present in the above-mentioned .csv files is in these files, in addition to parameters fit to the AOC models, eye-tracking information, and more. 
    allR: a structure with matrices containing statistics on behavioral performance for each individual subject. 
    rAvg: a structure similar to allR, but with each measure averaged across subjects. rAvg.valsByIndex is a guide to the dimensions of each matrix attached to allR and rAvg. 
For example, allR.PCs is a 5D matrix of the proportion correct achieved by each participant in each condition of the experiment. According to rAvg.valsByIndex, dimension 1 is the cue condition, where 0=distributed cue, 1=focal cue. Dimension 2 is the target side (the side the subject was required to judge on each trial), where 1=left, and 2=right. Dimension 3 is "congruent", which specifices whether the two words on each side came from the same category (congruent) or came from different categories (incongruent);  0=incongruent, and 1=congruent. Dimension 4 is "targPres", which specifies whether the word asked about on each trial was from the "living" category. "Living" words are categorized as "target present", for the purpose of signal detection theory analysis, and "non-living" words are "target absent". Dimension 5, the last dimension, is for individaul subjects (N=15). Note that the first entry of each dimension is for all the data *not* dividing by that variable. For example, to pull out all individual's proportion corrects on the left and right sides in the focal cue condition, you would execute this command: pcs = squeeze(allR.PCs(3,2:3,1,1,:)); That gives you a 2x15 matrix of p(correct) values, with one column for each subject and the first row for the left side of fixation and the second row for the rigth side of fixation. 

AllSubjBehavior_WordFreq.mat: 
    This file contains two variables, allR and rAvg, which are very much like the variables in AllSubjBehavior, except they are further sub-divided by the lexical frequency bin and length of each word. As shown in rAvg.valsByIndex, each matrix now has dimensions targLength and distLength, which are for the length of the target word (the one the subject had to judge) and the 'distractor' word (the word not asked about on that trial). There are also dimensions for targFreqBin  and distFreqBin (1=low, 2=high). 

AllSubjectVoxelResponseTable.csv:
This is a table of individual voxel data from the localizer scans and main experiment. Only voxels from the ROIs are included (for both hemispheres, V1, V2, V3, V4, VO, LO, VWFA-1, and VWFA-2). The voxel responses, expressed as percent signal change, were estimated by glmDenoise. The columns of this table are: 
        subject: integers 1-15, indicating subject numbers
        voxelIndex: linear index of each voxel within the subject's functional volume
        region: character string of the name of the ROI that this voxel belongs to. 
        hemisphere: 'Left' or 'Right'
        X, Y and Z: coordinates of this voxel in the functional volume, which aligned to ACPC space. 
        resp_wordL: response of this voxel to single words on the left of fixation in the localizer scans, as estimated by glmDenoise. Units = percent signal change. 
        resp_wordR: response of this voxel to single words on the right of fixation in the localizer scans.
        resp_scramL: response of this voxel to phase-scrambled word images on the left. 
        resp_scramR: response of this voxel to phase-scrambled word images on the right. 
        resp_focalCueLeft: response of this voxel in the focal cue left condition of the main experiment. 
        resp_focalCueRight: response of this voxel in the focal cue right condition of the main experiment.
        resp_distributedCue: response of this voxel in the distributed cue condition of the main experiment.
        xvalRSqr_locz: glmDenoises's statistic of cross-validated r-squared from the localizer scans 
        xvalRSqr_main: glmDenoises's statistic of cross-validated r-squared from the main experiment scans
        snrAfter_locz: glmDenoises's statistic of signal-to-noise ratio, after denoising, from the localizer scans.
        snrAfter_main: glmDenoises's statistic of signal-to-noise ratio, after denoising, from the main experiment scans.


AllSubjectVoxelResponseTable_WordFreq.csv: 
This is a table organized much like the table in AllSubjectVoxelResponseTable.csv, except the responses from the main experiment are divided by the lexical frequencies of the left and right words, which are each categorized as "high" or "low". This division is present for each cue condition: for example, for the distributed cue condition, each voxel has 4 responses: 
    resp_distrib_leftLow_rightLow, resp_distrib_leftLow_rightHigh, resp_distrib_leftHigh_rightLow, and resp_distrib_leftHigh_rightHigh. 


(3) MakeFiguresAndStats: this folder contains one Matlab script, AllFigureScript, and 21 separate Matlab function files, each of which makes a figure or prints some statistics. AllFigureScript calls all 21 of those functions in sequence. Figures are saved as .eps files in the Figures folder, and statistics are printed to .txt files in the Stats folder. 


(4) Utilities: this folder contains other Matlab files needed to run the analyses and create the figures. All of these functions were written by Alex White, except for: boyntonBootstrap (Geoffrey Boynton); exportfig (Ben Hinkle, available on the Mathworks FileExchange at https://www.mathworks.com/matlabcentral/fileexchange/727-exportfig); hist2D (David Bean; http://www.davidbdean.com/2006/08/17/how-to-plot-a-2d-histogram-using-matlab/).  


When you run the code, it creates additional folders: Figures, Results, and Stats. 

Finally, outside of any of the folders, is a getPaths.m function, which is used by many of the other functions to find the correct folders to read from and write to. 

This code was developed in Matlab2016b, using the Matlab Statistics toolbox, on an Apple Macbook Pro. 

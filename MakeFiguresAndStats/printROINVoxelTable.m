% Print information for Table 1 in White, Palmer, Boynton & Yeatman, PNAS
% 2019. 
% This table contains the mean numbers of voxels in each ROI. 
%
% Alex L. White
% University of Washington, 2019
% 
function printROINVoxelTable()

%% set paths
%add whole analysis code directory to the path
analysisDir = fileparts(fileparts(which(mfilename)));
addpath(genpath(analysisDir));

p = getPaths(); 

%% load data
resFile = fullfile(p.results,'AllSubjMainExptResponses.mat');
load(resFile);

%% open stats file 
sf = fopen(fullfile(p.stats, 'NVoxelTable.txt'),'w'); 

%% extract data 
areaIs = [2 3 4 5 6 1 7 8]; %rearrange visual areas slightly 
areaNames = allR.valsByIndex.brainArea(areaIs);

Ns = squeeze(allR.numVoxels(areaIs, :, :));

%% print mean number of voxels
fprintf(sf,'Number of voxels in each ROI\n'); 

fprintf(sf,'\tLH\tRH'); 
for ai=1:length(areaIs)
    fprintf(sf,'\n%s\t', areaNames{ai}); 
    for hi=1:2
        theseNs = squeeze(Ns(ai,hi,:)); 
        theseNs = theseNs(theseNs>0);
        mn = round(mean(theseNs)); 
        sem = round(standardError(theseNs'));        
        fprintf(sf,'%i (%i)\t', mn, sem);
    end
end

    
%% print N subjects with each ROI
fprintf(sf,'\n\nNumber of subjects with each ROI\n'); 

fprintf(sf,'\tLH\tRH'); 
for ai=1:length(areaIs)
    fprintf(sf,'\n%s\t', areaNames{ai}); 
    for hi=1:2
        theseNs = squeeze(Ns(ai,hi,:)); 
        nn = sum(theseNs>0); 
        fprintf(sf,'%i\t', nn);
    end
end

   

%Figure S5 in White, Palmer, Boynton & Yeatman, 2019
%Correlations between spatial selectivity and selective attention effects on
%individual voxels.
%
% by Alex L. White, University of Washington 2019
function FigS5_VoxelSpatialVsAttnSelectivity()

%% set paths
%add whole analysis code directory to the path
analysisDir = fileparts(fileparts(which(mfilename)));
addpath(genpath(analysisDir));

p = getPaths(); 

%% load big table
tableFileName = fullfile(p.data,'AllSubjectVoxelResponseTable.csv');
T = readtable(tableFileName);

%% open a file to print stats
statsFile = fullfile(p.stats,'StatsS5_SpatialVsAttnSelectivityStats.txt');
diary(statsFile);
statsF = fopen(statsFile,'w');

fprintf(1,'Analyais of relations between voxel spatial selectivity and selective attention effects\n');

%% extract data and set options
%regions to plot: all of retionotopic areas together, and VWFAs separately
regionsToPlot = {'Retinotopic','VWFA_1','VWFA_2'};
nRegions = length(regionsToPlot);
retinoI = find(strcmp(regionsToPlot,'Retinotopic'));
retinoRegions = {'V1','V2','V3','V4','LO','VO'};

%two hemispheres in different columns
hemNames = unique(T.hemisphere);

%% setup plot
xlimReti = [-2 2];
xlimVWFA = [-1.25 1.25];

axticksRetino = -2:1:2;
axticksVWFA = -1:0.5:1;

%bins for 2D histograms
nBins = 51;
minBin = -2.5; maxBin = 2.5;
binCenters = linspace(minBin, maxBin, nBins);

%set colormap
ncolrs = nBins;
hues = ones(1,ncolrs)*0.9;

minSat = 0.001; maxSat = 0.8;
minVal = 1; maxVal = 0.0001;

sats = linspace(minSat,maxSat,ncolrs);
vals = linspace(minVal,maxVal,ncolrs);

colormapColrs = hsv2rgb([hues' sats' vals']);

fontSize = 10;

%% compute max voxel count across all regions, for setting color map later
maxNs = zeros(nRegions, 2);
for ri = 1:nRegions
    for hemi = 1:2
        thisHem = hemNames{hemi};
        
        if ri==retinoI
            %add all the retinotopic areas
            roiSet = false(size(T.region));
            for rri=1:length(retinoRegions)
                roiSet = roiSet | strcmp(T.region, retinoRegions{rri});
            end
            roiSet = roiSet & strcmp(T.hemisphere, thisHem);
        else
            roiSet = strcmp(T.region,regionsToPlot{ri}) & strcmp(T.hemisphere, thisHem);
        end
        
        %x axis: hemifield selectivity during localizer:
        laterality = double(T.resp_wordL(roiSet,:) - T.resp_wordR(roiSet,:));
        
        %y axis: selective attention effect during main scan:
        selectiveAttn = double(T.resp_focalCueLeft(roiSet,:) - T.resp_focalCueRight(roiSet,:));
        
        [n,~,~] = hist2d(laterality,selectiveAttn,binCenters);
        maxNs(ri,hemi) =  max(n(:));
        
    end
end

%% plot and print stats

nCols=2;
nRows = nRegions;
figure;
for ri = 1:nRegions
    for hemi = 1:2
        subplot(nRows,nCols,(ri-1)*2+hemi); hold on;
        thisHem = hemNames{hemi};
        hemLabel = thisHem;
        
        if ri==retinoI
            %add all the retinotopic areas
            roiSet = false(size(T.region));
            for rri=1:length(retinoRegions)
                roiSet = roiSet | strcmp(T.region, retinoRegions{rri});
            end
            roiSet = roiSet & strcmp(T.hemisphere, thisHem);
            maxN = max(maxNs(ri,:));
            
            thisRegion = 'Retinotopic';
            xlims = xlimReti;
            axticks = axticksRetino;
        else %other regions separately
            roiSet = strcmp(T.region,regionsToPlot{ri}) & strcmp(T.hemisphere, thisHem);
            otherNs = maxNs(setdiff(1:nRegions,retinoI),:); %exclude all retinotopic
            maxN = max(otherNs(:));
            thisRegion = regionsToPlot{ri};
            xlims = xlimVWFA;
            axticks = axticksVWFA;
        end
        
        ylims = xlims;
        regionTitle = thisRegion;
        regionTitle(regionTitle == '_') = '-';
          
        %x axis: hemifield selectivity during localizer:
        laterality = double(T.resp_wordL(roiSet,:) - T.resp_wordR(roiSet,:));
        
        %y axis: selective attention effect during main scan:
        selectiveAttn = double(T.resp_focalCueLeft(roiSet,:) - T.resp_focalCueRight(roiSet,:));
        
        %correlation
        [corrRho, corrP] = corr(laterality, selectiveAttn);
        
        %also run corr for each subject
        subject = T.subject(roiSet,:);
        uSubjs = unique(subject);
        nSubjs = length(uSubjs);
        
        indivRhos = NaN(nSubjs,1);
        indivPs = NaN(nSubjs,1);
        for si=1:nSubjs
            subjSet = subject == uSubjs(si);
            if sum(subjSet)>0
                [indivRhos(si), indivPs(si)] = corr(laterality(subjSet), selectiveAttn(subjSet));
            end
        end
        
        meanRho = nanmean(indivRhos);          
        rhoCI = boyntonBootstrap(@mean, indivRhos(~isnan(indivRhos)),5000,95);
        
        fprintf(1,'\n--------------------------------------------------------------');
        fprintf(1,'\n--------------------------------------------------------------');
        fprintf(1,'\nANALYSIS OF %s %s', hemLabel, regionTitle);
        fprintf(1,'\n--------------------------------------------------------------');
        fprintf(1,'\n--------------------------------------------------------------\n');
         
        fprintf(1,'\nPearson linear correlation between laterality and selective attention:\n');
        fprintf(1,'Union of all subjects'' voxels: rho=%.4f, p=%.4f\n',corrRho, corrP);
        fprintf(1,'\nCorrelations on each individual subject:\n');
        fprintf(1,'mean rho = %.3f,  95%% CI = [%.3f %.3f]\n', meanRho, rhoCI(1), rhoCI(2));
        
        %make a new table to model this region
        d = table(subject, laterality, selectiveAttn);
        
        %1. Simplest model with no random effects: 
        eqtn1a = 'selectiveAttn ~ laterality';
        nosubjectModel = fitlme(d,eqtn1a);
        
        %2. More complex model with has random intercepts for subjects 
        eqtn2 = 'selectiveAttn ~ laterality + (1 | subject)';
        lme2 = fitlme(d,eqtn2);
        
        %3. Most complex model with random intercepts and slope (effect of laterality) for subjects
        eqtn3 = 'selectiveAttn ~ laterality + (laterality | subject)';
        lme3 = fitlme(d,eqtn3);
        
        %model comparisoin
        fprintf(1,'\nComparison of models with and without random slopes by subject:\n');
        comp2 = compare(lme2,lme3)
        compT2 = dataset2table(comp2);
        compP2 = compT2.pValue(end);
        
        %if the bigger model with random slopes by subject is significantly better, then use it
        if compP2<0.05
            subjectModel  = lme3;
        else
            subjectModel  = lme2;
        end
        
        %compare the better of the two big models with the simple model with no random effects at all
        fprintf(1,'\nComparison of the simplest model without random effects against the wining model with subject effects:\n');
        comp = compare(nosubjectModel,subjectModel)
        compT = dataset2table(comp);
        compP = compT.pValue;
        
        if compP<0.05
            bestModel  = subjectModel; 
            if compP2<0.05
                fprintf(1,'\nFor %s %s, we used the LME model with random intercepts AND slopes by subject\n', hemLabel, regionTitle);
            else
                fprintf(1,'\nFor %s %s, we used the LME model with random intercepts (but not slopes) by subject\n', hemLabel, regionTitle);
            end
        else
            bestModel = nosubjectModel;
            fprintf(1,'\nFor %s %s, we used the simplest LME model with no random effects\n', hemLabel, regionTitle);
        end
        
        display(bestModel);
        
        %compute prediction from the best-fitting model:
        [~,~,lmeStats] = fixedEffects(bestModel);
        lmeIntercept = lmeStats.Estimate(1);
        lmeSlope = lmeStats.Estimate(2);     
        
        predXs = [min(laterality) max(laterality)];
        predXs(predXs>max(xlims)) = max(xlims);
        predXs(predXs<min(xlims)) = min(xlims);
        
        lmePredYs = lmeIntercept + predXs*lmeSlope;        
       
        %% plot
        %2D histogram
        [n,x,y] = hist2d(laterality,selectiveAttn,binCenters);
        
        %set color map
        colormap(colormapColrs);
        colorClip = [0 maxN];
        
        %image the 2D histogram
        imagesc(x(1,:),y(:,1),n,colorClip);
        
        if hemi==2
            cbh=colorbar;
            cbh.Label.String = 'Num. voxels';
        end
        
        plot([0 0],ylims,'k--');
        plot(xlims,[0 0],'k--');
        
        %plot overall fit
        plot(predXs,lmePredYs,'r-');
        
        %print the correlation coefficient and 95% CI 
        textx = xlims(1)+0.085*diff(xlims);
        texty = ylims(1)+0.1*diff(ylims);
        text(textx,texty, sprintf('rho=%.2f, [%.2f %.2f]', meanRho, rhoCI(1), rhoCI(2)));
        
        axis square;
        
        xlim(xlims); ylim(ylims);
        set(gca,'XTick',axticks,'YTick',axticks);
        
        if ri==nRows
            xlabel('left word - right word');
        end
        if hemi==1
            ylabel('attnd L - attnd R');
        end
        
        title([hemLabel ' ' regionTitle]);
    end
end

diary off;

figName = 'FigS5_VoxelSpatialVsAttnSelectivity.eps';

set(gcf,'units','centimeters','pos',[5 5 17 5.9*nRegions]);
exportfig(gcf,fullfile(p.figures,figName),'Format','eps','bounds','loose','color','rgb','LockAxes',0,'FontMode','fixed','FontSize',fontSize);
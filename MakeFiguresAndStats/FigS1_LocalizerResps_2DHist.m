%Figure S1 in White, Palmer, Boynton & Yeatman, 2019
%2D histograms of individual voxel responses to single words on the left vs
%singe words on the right (from localizer scan). 
%
% by Alex L. White, University of Washington 2019
function FigS1_LocalizerResps_2DHist()

%% set paths
%add whole analysis code directory to the path
analysisDir = fileparts(fileparts(which(mfilename)));
addpath(genpath(analysisDir));

p = getPaths(); 
%% load big table

tableFileName = fullfile(p.data,'AllSubjectVoxelResponseTable.mat');
load(tableFileName,'T');


%% set which data to plot 
%regions to plot: all of retionotopic areas together, and VWFAs separately
regionsToPlot = {'Retinotopic','VWFA_1','VWFA_2'};
nRegions = length(regionsToPlot);
retinoI = find(strcmp(regionsToPlot,'Retinotopic'));
retinoRegions = {'V1','V2','V3','V4','LO','VO'};

%and hemispheres in separate columns
hemTitles = unique(T.hemisphere);

uSubjs = unique(T.subject);
nSubj = numel(uSubjs);

%whether to plot dots over each subject's median responses
plotSubjMedians = true;

%% plot parameters
fontSize = 10;

axlim = [-1 2.5];
axticks = -1:1:2;

nBins = 51;
binCenters = linspace(axlim(1), axlim(2), nBins);

%set colormap
ncolrs = nBins;
hues = ones(1,ncolrs)*0.9;

minSat = 0.001; maxSat = 0.8;
minVal = 1; maxVal = 0.0001;

sats = linspace(minSat,maxSat,ncolrs);
vals = linspace(minVal,maxVal,ncolrs);

colormapColrs = hsv2rgb([hues' sats' vals']);

axLineColr = 'k';
subjMedianColr='r';

nCols = 2;
nRows = nRegions;

%% first, find the max voxel count in each region, so we can set the colormap appropriately 
maxNs = zeros(nRegions, 2);
for ri = 1:nRegions
    for hemi = 1:2
        thisHem = hemTitles{hemi};
        
        if ri==retinoI
            %add all the retinotopic areas
            subst = false(size(T.region));
            for rri=1:length(retinoRegions)
                subst = subst | strcmp(T.region,retinoRegions{rri});
            end
            subst = subst & strcmp(T.hemisphere, thisHem);
        else
            thisRegion = regionsToPlot{ri};
            subst = strcmp(T.region,thisRegion) & strcmp(T.hemisphere, thisHem);
        end
        
        xs = T.resp_wordL(subst,:);
        ys = T.resp_wordR(subst,:);
        [n,~,~] = hist2d(xs,ys,binCenters);
        maxNs(ri,hemi) =  max(n(:));
    end
end

%% then plot 
figure;
for ri = 1:nRegions
    for hemi = 1:2
        thisHem = hemTitles{hemi};
        
        subplot(nRows,nCols,(ri-1)*2+hemi); hold on;
        
        %gather all the retinotopic areas together 
        if ri==retinoI
            subst = false(size(T.region));
            for rri=1:length(retinoRegions)
                subst = subst |  strcmp(T.region,retinoRegions{rri});
            end
            subst = subst & strcmp(T.hemisphere, thisHem);
            regionTitle = 'Retinotopic';
            
            maxN = max(maxNs(ri,:));
        %other areas separately     
        else
            thisRegion = regionsToPlot{ri};
            subst = strcmp(T.region,thisRegion) & strcmp(T.hemisphere, thisHem);
            regionTitle = thisRegion;
            
            otherNs = maxNs(setdiff(1:nRegions,retinoI),:); %exclude all retinotopic
            maxN = max(otherNs(:));
            
        end
        regionTitle(regionTitle == '_') = '-';
        
        subjects = T.subject(subst);
        
        %pull out x and y data: reponses to single words on the left and right 
        xs = T.resp_wordL(subst,:);
        ys = T.resp_wordR(subst,:);
        
        %compute the 2D histograpm
        [n,x,y] = hist2d(xs,ys,binCenters);

        %clip the color at max for this ROI
        colorClip = [0 maxN];
        colormap(colormapColrs);
        
        %image the 2D histogram: 
        imagesc(x(1,:),y(:,1),n,colorClip);
        
        if hemi==2
            cbh=colorbar;
            cbh.Label.String = 'Num. voxels';
        end
        
        %put red dots at individual subject medians
        if plotSubjMedians
            for subji=1:nSubj
                subjectSet = subji == subjects;
                plot(median(xs(subjectSet)),median(ys(subjectSet)),'.','MarkerSize',6, 'Color', subjMedianColr);
            end
        end
        
        %plot 0 axes
        plot([0 0],axlim,'--','Color',axLineColr);
        plot(axlim,[0 0],'--','Color',axLineColr);
        
        %and the identity line
        plot(axlim, axlim, '-','Color',axLineColr);
        
        xlim(axlim); ylim(axlim);
        axis square;
        set(gca,'XTick',axticks,'YTick',axticks);
        
        if ri==nRegions
            xlabel('Left word (psc)');
        end
        if hemi==1
            ylabel('Right word (psc)');
        end        
        title([hemTitles{hemi} ' ' regionTitle]);
    end
end

%save figure;
figName = sprintf('FigS1_Localizer2DHist.eps');
set(gcf,'units','centimeters','pos',[5 5 17 5.9*nRegions]);
exportfig(gcf,fullfile(p.figures,figName),'Format','eps','bounds','loose','color','rgb','LockAxes',0,'FontMode','fixed','FontSize',fontSize);


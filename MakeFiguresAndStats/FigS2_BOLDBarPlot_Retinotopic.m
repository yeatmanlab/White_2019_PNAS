%% Figure S2 for White, Palmer, Boynton & Yeatman, PNAS 2019
% Makes a bar plot of mean BOLD responses in the main attention experiment, 
% in individual retinotopic areas only. 
% In the second and third subplots, it shows 
% the mean selective and divided attention effects on mean BOLD responses. 
% 
% by Alex L. White at the University of Washington, 2019
%
function FigS2_BOLDBarPlot_Retinotopic()

%% set paths
%add whole analysis code directory to the path
analysisDir = fileparts(fileparts(which(mfilename)));
addpath(genpath(analysisDir));

p = getPaths();

%% load data
resFile = fullfile(p.results,'AllSubjMainExptResponses.mat'); 
load(resFile); 

%% extract data 

condIs = 1:3;
condLabs = allR.valsByIndex.condition(condIs);
nConds = numel(condLabs);

hemIs = 1:2;
hemLabs = allR.valsByIndex.hemisphere(hemIs);
nHems = length(hemIs);

%Average over Retinotopic cortex
areaNames = {'V1','V2','V3','V4','LO','VO'};
[~,areasToPlot] = intersect(allR.valsByIndex.brainArea,areaNames,'stable');

%reorder these
areasToPlot = areasToPlot([2 3 4 5 6 1]);

betas = allR.meanBetas;
nAreas = length(areasToPlot);

%label the brain areas and take out hyphens
areaLabs = allR.valsByIndex.brainArea(areasToPlot);
for ai=1:nAreas
    al = areaLabs{ai};
    al(al=='_') = '-';
    areaLabs{ai} = al;
end

nSubj = size(betas,ndims(betas));

%% compute means and SEms
allbs = squeeze(betas(areasToPlot,hemIs,condIs,:));
meanBs = nanmean(allbs,ndims(allbs));
semBs = standardError(allbs,ndims(allbs));

%% ATTENTION EFFECTS
nComps = 2;
ds       = NaN(nAreas,nHems,nComps,nSubj);
bootCIs  = NaN(nAreas, nHems, nComps,2);
bootSigs = NaN(nAreas, nHems, nComps);
bootPs   = NaN(nAreas, nHems, nComps);

focalLeftI  = find(strcmp(condLabs,'Focal cue left')); 
focalRightI = find(strcmp(condLabs,'Focal cue right'));
distribI    = find(strcmp(condLabs,'Distributed cue')); 

effectTitles{1} = 'Focal contralateral - focal ipsilateral';
effectTitles{2} = 'Focal contralateral - distributed';

CIBiasCorrect = false;

for hi=1:2 
    if strcmp(hemLabs{hi},'Left')
        cueComps = [focalRightI focalLeftI; focalRightI distribI]; %contra - ipsi cued; contra cued - distributed;
    elseif strcmp(hemLabs{hi},'Right')
        cueComps = [focalLeftI focalRightI; focalLeftI distribI]; %contra - ipsi cued; contra cued - distributed;
    end
    
    for ei=1:nComps
        ds(:,hi,ei,:) = squeeze(allbs(:,hi,cueComps(ei,1),:) - allbs(:,hi,cueComps(ei,2),:));
    end
    for ai=1:nAreas
        for ei=1:nComps
            [diffCI,~,bootDist] = boyntonBootstrap(@nanmean,squeeze(ds(ai,hi,ei,:))',5000,95,CIBiasCorrect);
            bootSigs(ai,hi,ei) = all(diffCI>0) || all(diffCI<0);
            bootCIs(ai,hi,ei,:) = diffCI;
            bootPs(ai,hi,ei) = getBootPs(bootDist',0);
        end
    end
end
meanDs = nanmean(ds,ndims(ds));
semDs  = standardError(ds,ndims(ds));

%% set up plot
barWidth = 0.08;
edgeLineWidth = 1;
effectBarWid = 0.09;
semLineWid = 1;
axLineWidth = 1;
fontSize = 9;

xjitt = [-0.29 -0.18 -0.07;  0.07 0.18 0.29];

xlims = [0.5 nAreas+0.5];
ylims = [-0.25 0.5];
legendLoc = 'NorthWest';

diffYLims = [-0.1 0.15];
effectXJitt = [-0.11 0.11];

cueHues = [0.8 0.6 0.4];
cueSats = [0.9 0.9 0.9];
cueVals = [0.7 0.7 0.7];
cueColrs = hsv2rgb([cueHues' cueSats' cueVals']);

fillColors = zeros(2,3,3);
fillColors(1,:,:) = cueColrs;
fillColors(2,:,:) = ones(3,3);

edgeColors = zeros(2,3,3);
edgeColors(1,:,:) = cueColrs;
edgeColors(2,:,:) = cueColrs;

selectiveColor = hsv2rgb([0.5 0.9 0.6]);
dividedColor   = hsv2rgb([0.7 0.9 0.85]);

diffErrorBarType = 'sem'; %or CI

figHeightCm = 14;
MRIWidCm = 17.8;

%% variable subplot sizes and positions
leftMargin = 0.1;
rightMargin = 0.04;
topMargin = 0.025;
bottomMargin = 0.07;

subplotVerticalSpacing = 0.085;
subplotHorizontalSpacing = 0.05;

nRows = 3; nCols = 1;

subplotRelativeHeights = [2 1 1];

availableHeight =  (1-topMargin-bottomMargin-subplotVerticalSpacing*(nRows-1));
subplotHeights = availableHeight*subplotRelativeHeights/sum(subplotRelativeHeights);
subplotWidths = ones(1,nCols)*(1-leftMargin-rightMargin);

subplotPositions = zeros(nRows,nCols,4);

for ri=1:nRows
    for ci=1:nCols
        leftPos = leftMargin + (ci-1)*subplotHorizontalSpacing + sum(subplotWidths(1:(ci-1)));
        
        plotsBelow = (ri+1):nRows;
        nVertSpacing = length(plotsBelow);
        if nVertSpacing<0, nVertSpacing = 0; end
        bottomPos = bottomMargin + nVertSpacing*subplotVerticalSpacing + sum(subplotHeights(plotsBelow));
        
        subplotPositions(ri,ci,:) = [leftPos bottomPos subplotWidths(ci) subplotHeights(ri)];
    end
end


figure; hold on;
%% 1. plot mean BOLD responses in 3 attention conditions

rowI = 1; colI = 1;
subplot('position',squeeze(subplotPositions(rowI,colI,:)));
hold on;

condHs = zeros(1,nConds);

plot(xlims,[0 0],'k-');

for ai = 1:nAreas
    for hi=1:nHems
        xs = ai+xjitt(hi,:);
        for ci=1:nConds
            bx = xs(ci)+[-0.5 0.5]*barWidth;
            by = [0 meanBs(ai,hi,ci)];
            
            vertx=[bx; bx];
            verty=[by fliplr(by)];
            
            %bar
            ahand = fill(vertx(:), verty(:), squeeze(fillColors(hi,ci,:))', 'EdgeColor', squeeze(edgeColors(hi,ci,:))', 'LineWidth', edgeLineWidth);
            if hi==1, condHs(ci) = ahand; end
            
            %do it again in reverse order to avoid diagonal white line
            fill(flipud(vertx(:)), flipud(verty(:)), squeeze(fillColors(hi,ci,:))', 'EdgeColor', squeeze(edgeColors(hi,ci,:))', 'LineWidth', edgeLineWidth);
            
            %error bar
            plot(xs([ci ci]),meanBs(ai,hi,ci)+[-1 1]*semBs(ai,hi,ci),'k-','LineWidth',semLineWid);
        end
    end
end

xlim(xlims);
ylim(ylims);

yticks = ylims(1):0.25:ylims(2);

set(gca,'XTick',1:nAreas,'YTick',yticks,'LineWidth',axLineWidth);
set(gca,'XTickLabels',areaLabs);
ylabel('percent signal change (psc)');
legend(condHs,condLabs,'Location',legendLoc);  

%% 2 & 3. Selective and Divided attention effects
for ei=1:2
    if ei==1
        effectFillColors  = [selectiveColor; 1 1 1];
        effectEdgeColors = [selectiveColor; selectiveColor];
        
    else
        effectFillColors  = [dividedColor; 1 1 1];
        effectEdgeColors = [dividedColor; dividedColor];
        
    end
    
    rowI = ei+1; colI = 1;
    subplot('position',squeeze(subplotPositions(rowI,colI,:)));
    hold on;
    
    plot(xlims,[0 0],'k-');
    for ai = 1:nAreas
        effectHands = zeros(1,nHems);
        
        for hi = 1:nHems
            xs = ai+effectXJitt(hi);
            
            bx = xs+[-0.5 0.5]*effectBarWid;
            by = [0 meanDs(ai,hi,ei)];
            
            vertx=[bx; bx];
            verty=[by fliplr(by)];
            
            effectHands(hi) = fill(vertx(:), verty(:), effectFillColors(hi,:), 'EdgeColor', effectEdgeColors(hi,:), 'LineWidth', edgeLineWidth);
            %do it again in reverse order to avoid diagonal white line
            fill(flipud(vertx(:)), flipud(verty(:)), effectFillColors(hi,:), 'EdgeColor', effectEdgeColors(hi,:), 'LineWidth', edgeLineWidth);
            
            %plot error bar
            if strcmp(diffErrorBarType,'sem')
                plot(xs([1 1]),meanDs(ai,hi,ei)+[-1 1]*semDs(ai,hi,ei),'k-','LineWidth',semLineWid);
                
                %star if significant by bootstrapping
                thisP =  bootPs(ai,hi,ei);
                if thisP<0.05
                    starY = diffYLims(1)+0.9*diff(diffYLims);
                    if thisP<0.001
                        startTxt = '***';
                    elseif thisP<0.01
                        startTxt = '**';
                    elseif thisP<0.05
                        startTxt = '*';
                    end
                    text(xs(1),starY,startTxt,'HorizontalAlignment','center');
                end
            else %95% CI
                plot(xs([1 1]),squeeze(bootCIs(ai,hi,ei,:)),'k-','LineWidth',semLineWid);
            end
        end
    end
    
    xlim(xlims);
    ylim(diffYLims);
    set(gca,'XTick',1:nAreas,'YTick',diffYLims(1):0.05:diffYLims(2),'LineWidth',axLineWidth);
    yticklabels = get(gca,'YTickLabel');
    yticklabels(2:2:end) = {''};
    set(gca,'YTickLabel',yticklabels);
    
    if ei==2, set(gca,'XTickLabels',areaLabs);
    else, set(gca,'XTickLabels',{}); end
    ylabel('\Delta psc');
    title(effectTitles{ei});
    legend(effectHands,hemLabs);
end

%save figure;
set(gcf,'units','centimeters','pos',[8 8 MRIWidCm figHeightCm],'color','w');
figTitle = fullfile(p.figures,'FigS2_BOLDAttnEffects_Retinotopic.eps');
exportfig(gcf,figTitle,'Format','eps','bounds','loose','color','rgb','LockAxes',0,'FontMode','fixed','FontSize',fontSize);


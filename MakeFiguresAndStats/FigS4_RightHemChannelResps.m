%% Figure S4 for White, Palmer, Boynton & Yeatman, PNAS 2019
% Makes a bar plot of mean channel responses in *right* hemisphere VWFAs.
% In the second and third subplots, it plots  the mean selective and
% divided attention effects on channel responses.
%
% by Alex L. White at the University of Washington, 2019
%
function FigS4_RightHemChannelResps()

%% set paths
%add whole analysis code directory to the path
analysisDir = fileparts(fileparts(which(mfilename)));
addpath(genpath(analysisDir));

p = getPaths(); 

%% load data 
resFile = fullfile(p.results,'AllSubjChannelResponses.mat');
load(resFile);

%% pull out data
areasToPlot = [find(strcmp(allR.valsByIndex.brainArea,'VWFA_1')) find(strcmp(allR.valsByIndex.brainArea,'VWFA_2'))];
nAreas = length(areasToPlot);
brainAreas = allR.valsByIndex.brainArea(areasToPlot);

hemI = find(strcmp(allR.valsByIndex.hemisphere,'Right'));
nHems = length(hemI);
hemLab = allR.valsByIndex.hemisphere{hemI};

channIs = 1:2;
channelLabels = allR.valsByIndex.channel(channIs);
nChann = length(channIs);
for ki=1:nChann
    channelLabels{ki} = sprintf('%s channel', channelLabels{ki});
end

condsToPlot = 1:3;
condLabs = allR.valsByIndex.condition(condsToPlot);
nCond = numel(condLabs);

areaLabs = cell(1,nAreas);
%take out hyphens
for ai=1:nAreas
    al = brainAreas{ai};
    al(al=='_') = '-';
    areaLabs{ai} = al;
end

rs = squeeze(allR.twoChannelResponses(areasToPlot,hemI,channIs,condsToPlot,:));
meanResps = nanmean(rs,ndims(rs));
semResps = standardError(rs,ndims(rs));

nSubj = size(rs,ndims(rs));

% ATTENTION EFFECTS
nComps = 2;
ds = NaN(nAreas,nChann,nComps,nSubj);
bootCIs = NaN(nAreas, nChann, nComps,2);
bootSigs = NaN(nAreas, nChann, nComps);
bootPs = NaN(nAreas, nChann, nComps);

focalLeftI  = find(strcmp(condLabs,'Focal cue left'));
focalRightI = find(strcmp(condLabs,'Focal cue right'));
distribI    = find(strcmp(condLabs,'Distributed cue'));

CIBiasCorrect = false;

for ki=1:nChann %left channel, right channel
    if strcmp(channelLabels{ki},'Left channel')
        cueComps = [focalLeftI focalRightI; focalLeftI distribI]; %focal cued - uncued; focal cued  - distributed
    elseif strcmp(channelLabels{ki},'Right channel')
        cueComps = [focalRightI focalLeftI; focalRightI distribI]; %focal cued - uncued; focal cued  - distributed
    end
    
    for ei=1:nComps
        ds(:,ki,ei,:) = squeeze(rs(:,ki,cueComps(ei,1),:) - rs(:,ki,cueComps(ei,2),:));
    end
    
    for ai=1:nAreas
        for ei=1:nComps
            [diffCI,~,bootDist] = boyntonBootstrap(@nanmean,squeeze(ds(ai,ki,ei,:))',5000,95,CIBiasCorrect);
            bootSigs(ai,ki,ei) = all(diffCI>0) || all(diffCI<0);
            bootCIs(ai,ki,ei,:) = diffCI;
            bootPs(ai,ki,ei) = getBootPs(bootDist',0);
        end
    end
end

effectTitles{1} = 'Focal cued - uncued';
effectTitles{2} = 'Focal cued - distributed cue';

meanDs = nanmean(ds,ndims(ds));
semDs  = standardError(ds,ndims(ds));

%% set up plot
barWidth = 0.08;
edgeLineWidth = 1;
effectBarWid = 0.09;
semLineWid = 1;
axLineWidth = 1;
fontSize = 10;

xlims = [0.5 nAreas+0.5];
ylims = [0 1.75];
ytickd = 0.25;
diffYLims = [-0.2 0.4];
diffytickd = 0.2;

xjitt = [-0.31 -0.2 -0.09;  0.09 0.2 0.31];
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

figHeightCm = 14;
MRIWidCm = 8.9;

legendLoc = 'NorthWest';
diffErrorBarType = 'sem'; %or '95CI';

%% variable subplot sizes and positions
leftMargin = 0.175;
rightMargin = 0.035;
topMargin = 0.05;
bottomMargin = 0.07;

subplotVerticalSpacing = 0.085;
subplotHorizontalSpacing = 0.045;

nRows = 3; nCols = nHems;
subplotRelativeHeights = [2 1 1];

availableHeight =  (1-topMargin-bottomMargin-subplotVerticalSpacing*(nRows-1));
subplotHeights = availableHeight*subplotRelativeHeights/sum(subplotRelativeHeights);

subplotRelativeWidths = ones(1,nHems);
availableWidth = (1-leftMargin - rightMargin - subplotHorizontalSpacing*(nCols-1));
subplotWidths = availableWidth*subplotRelativeWidths/sum(subplotRelativeWidths);

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


%% Panel A. Mean Channel responses in 3 attention conditions
figure;
rowI = 1; colI = 1;
subplot('position',squeeze(subplotPositions(rowI,colI,:)));
hold on;

condHs = zeros(1,nCond);
plot(xlims,[0 0],'k-');

for ai = 1:nAreas
    for ki=1:nChann
        xs = ai+xjitt(ki,:);
        for ci=1:nCond
            bx = xs(ci)+[-0.5 0.5]*barWidth;
            by = [0 meanResps(ai,ki,ci)];
            
            vertx=[bx; bx];
            verty=[by fliplr(by)];
            
            %bar
            ahand = fill(vertx(:), verty(:), squeeze(fillColors(ki,ci,:))', 'EdgeColor', squeeze(edgeColors(ki,ci,:))', 'LineWidth', edgeLineWidth);
            %do it again in reverse order to avoid diagonal white line
            fill(flipud(vertx(:)), flipud(verty(:)), squeeze(fillColors(ki,ci,:))', 'EdgeColor', squeeze(edgeColors(ki,ci,:))', 'LineWidth', edgeLineWidth);

            if ki==1, condHs(ci) = ahand; end

            %error bar
            plot(xs([ci ci]),meanResps(ai,ki,ci)+[-1 1]*semResps(ai,ki,ci),'k-','LineWidth',semLineWid);
        end
    end
end

xlim(xlims);
ylim(ylims);

yticks = ylims(1):ytickd:ylims(2);

set(gca,'XTick',1:nAreas,'XTickLabels',areaLabs,'YTick',yticks,'LineWidth',axLineWidth);
ylabel('response');

legend(condHs,condLabs,'Location',legendLoc); 


%% Panels B & C: Selective and Divided attention effects
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
        for ki = 1:nChann
            xs = ai+effectXJitt(ki);
            
            bx = xs+[-0.5 0.5]*effectBarWid;
            by = [0 meanDs(ai,ki,ei)];
            
            vertx=[bx; bx];
            verty=[by fliplr(by)];
            
            fill(vertx(:), verty(:), effectFillColors(ki,:), 'EdgeColor', effectEdgeColors(ki,:), 'LineWidth', edgeLineWidth);
            %do it again in reverse order to avoid diagonal white line
            fill(flipud(vertx(:)), flipud(verty(:)), effectFillColors(ki,:), 'EdgeColor', effectEdgeColors(ki,:), 'LineWidth', edgeLineWidth);
            
            %plot error bar
            if strcmp(diffErrorBarType,'sem')
                plot(xs([1 1]),meanDs(ai,ki,ei)+[-1 1]*semDs(ai,ki,ei),'k-','LineWidth',semLineWid);
                
                %star if significant by bootstrapping
                thisP =  bootPs(ai,ki,ei);
                if thisP<0.05
                    starY = diffYLims(1)+0.95*diff(diffYLims);
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
                plot(xs([1 1]),squeeze(bootCIs(ai,ki,ei,:)),'k-','LineWidth',semLineWid);
            end
        end
    end
    
    xlim(xlims);
    ylim(diffYLims);
    set(gca,'XTick',1:nAreas,'YTick',diffYLims(1):diffytickd:diffYLims(2),'LineWidth',axLineWidth);
    
    if ei==2, set(gca,'XTickLabels',areaLabs);
    else, set(gca,'XTickLabels',{}); end
    ylabel('\Delta response');
    
    title(effectTitles{ei});
end
%save figure;
set(gcf,'units','centimeters','pos',[8 8 MRIWidCm figHeightCm],'color','w');

figTitle = fullfile(p.figures,sprintf('FigS4_ChannelResps_%sHem.eps',hemLab));
exportfig(gcf,figTitle,'Format','eps','bounds','loose','color','rgb','LockAxes',0,'FontMode','fixed','FontSize',fontSize);


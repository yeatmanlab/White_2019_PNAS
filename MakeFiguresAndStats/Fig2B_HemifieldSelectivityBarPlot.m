%% Figure 2B for White, Palmer, Boynton & Yeatman, PNAS 2019
% Makes a bar plot of mean BOLD responses in the localizer scan, to show 
% the hemifield selectivity. 
% Plots responses in the average of retinotopic areas, and in VWFA-1
% and VWFA-2, in both hemispheres. 
%
% Additionally, in a separate figure, it plots the mean lateralization
% indices (1-ipsi/contra). 
% 
% 
% by Alex L. White at the University of Washington, 2019
%
function Fig2B_HemifieldSelectivityBarPlot()

%% set paths
%add whole analysis code directory to the path
analysisDir = fileparts(fileparts(which(mfilename)));
addpath(genpath(analysisDir));

p = getPaths(); 

%% load data

resFile = fullfile(p.results,'AllSubjLocalizerResponses.mat'); 
load(resFile); 

%% extract data 

condsToPlot = {'left_word','right_word'};

nConds = numel(condsToPlot);
condIs = zeros(1,nConds);
condLabs = cell(1,nConds);
for ci=1:nConds
    condIs(ci) = find(strcmp(allR.valsByIndex.condition,condsToPlot{ci}));
    rl = condsToPlot{ci};
    rl(rl=='_') = ' ';
    condLabs{ci} = rl;
end

hemIs = 1:2;
hemLabs = allR.valsByIndex.hemisphere(hemIs);

%Average over Retinotopic cortex
areasToAverage = {'V1','V2','V3','V4','LO','VO'};
[~,areaToAverageIs] = intersect(allR.valsByIndex.brainArea,areasToAverage,'stable');
retinoBeta = nanmean(allR.meanBetas(areaToAverageIs,:,:,:),1);

%Append the mean responses in retinotopic areas to the matrix of betas
betas = cat(1,allR.meanBetas,retinoBeta);

%figure out which areas to plot
retinoI = size(betas,1);
brainAreas = cat(1,allR.valsByIndex.brainArea,'Retinotopic');
areasToPlot = [retinoI find(strcmp(brainAreas,'VWFA_1')) find(strcmp(brainAreas,'VWFA_2'))];
nAreas = length(areasToPlot);

%label the brain areas and take out hyphens
areaLabs = brainAreas(areasToPlot);
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


%% compute lateralization index: 1-(ipsilateral/contralateral)
allds = NaN(length(areasToPlot), length(hemIs), nSubj);
for hi=1:2
    % for LH, the ratio is right/left
    % for RH, the ratio is left/right
    allds(:,hi,:) = 1 - (allbs(:,hi,hi,:)./allbs(:,hi,3-hi,:));
end

mds = nanmean(allds,ndims(allds));
semDs = standardError(allds,ndims(allds));

%stats on those indices
dCIs = NaN(nAreas,2,2);
bootPs = NaN(nAreas,2);
for ai=1:nAreas
    for hi=1:2
        theseds = squeeze(allds(ai,hi,:));
        theseds = theseds(~isnan(theseds));
        [diffCI, ~,bootDist] = boyntonBootstrap(@nanmean,theseds,5000,95);
        dCIs(ai,hi,:) = diffCI;
        bootPs(ai,hi) = getBootPs(bootDist',0);
    end
end

%% plot parameters

xjitt = [-.25 -0.125; 0.125 0.25];
xlims = [0.5 nAreas+0.5];
ylims = [-0.2 1];
yticks = [0 0.5 1];

bwid = 0.09;
semLineWid = 1;

%bar colors: blue & green
edgeHSVs = [0.6 0.9 0.7;
            0.4 0.9 0.7];

edgeColrs = hsv2rgb(edgeHSVs);

fillColrs{1} = edgeColrs; %left hemisphere
fillColrs{2} = ones(size(edgeColrs)); %right hemisphere

errorBarColrs = edgeColrs/2;

fontSize = 9;
figWidth = 8; 
figHeight = 4;

%% plot mean BOLD responses
figure; hold on;

condHs = zeros(1,nConds);

plot(xlims,[0 0],'k-');

for ai=1:nAreas
    for hi = 1:2
        xs = ai+xjitt(hi,:);
        
        for ci=1:nConds
            %bar plot
            bx = xs(ci)+[-bwid/2 bwid/2];
            by = [0 meanBs(ai,hi,ci)];
            
            vertx=[bx; bx];
            verty=[by fliplr(by)];
            
            condHs(hi,ci) = fill(vertx(:), verty(:),fillColrs{hi}(ci,:),'EdgeColor',edgeColrs(ci,:),'LineWidth',1);
            %do it again in reverse order to avoid diagonal white line
            fill(flipud(vertx(:)), flipud(verty(:)),fillColrs{hi}(ci,:),'EdgeColor',edgeColrs(ci,:),'LineWidth',1);
          
            %error bars: +/- 1 SEM
            plot(xs([ci ci]),meanBs(ai,hi,ci)+[-1 1]*semBs(ai,hi,ci),'-','Color',errorBarColrs(ci,:),'LineWidth',semLineWid);
        end
    end
end

set(gca,'XTick',1:nAreas,'XTickLabels',areaLabs,'YTick',yticks);

xlim(xlims);
ylim(ylims); 

ylabel('% signal change');

legend(condHs(1,:),condLabs,'Location','NorthWest'); %legend boxoff;

%save figure;
set(gcf,'units','centimeters','pos',[8 8 figWidth figHeight],'color','w');
figTitle = fullfile(p.figures,'Fig2B_LocalizerWordResps.eps');
exportfig(gcf,figTitle,'Format','eps','bounds','loose','color','rgb','LockAxes',0,'FontMode','fixed','FontSize',fontSize);


%% Plot mean lateralization indices between stimulus left and right

figure; hold on;

diffJitt = [-0.1 0.1];

diffEdgeColr = hsv2rgb([0.5 0.9 0.6]);
diffFillColr = [diffEdgeColr; 1 1 1];

ylims = [0 1.5];
yLab = '1 - ipsi / contra';

plot(xlims,[0 0],'k-');

for ai=1:nAreas
    xs = ai+diffJitt;
    hemHs = NaN(1,2);
    for hi=1:2
        ds = squeeze(mds(ai,hi));
        %error bars: standard error
        es = ds+[-1 1]*squeeze(semDs(ai,hi));
        
        %bar plot
        bx = xs(hi)+[-bwid/2 bwid/2];
        by = [0 ds];
        
        vertx=[bx; bx];
        verty=[by fliplr(by)];
        
        hemHs(hi) = fill(vertx(:), verty(:),diffFillColr(hi,:),'EdgeColor',diffEdgeColr,'LineWidth',1);
        %do it again in reverse order to avoid diagonal white line
        fill(flipud(vertx(:)), flipud(verty(:)),diffFillColr(hi,:),'EdgeColor',diffEdgeColr,'LineWidth',1);
        
        plot(xs([hi hi]),es,'-','Color','k','LineWidth',semLineWid);
        
        thisP =  bootPs(ai,hi);
        if thisP<0.05
            starY = ylims(1)+0.9*diff(ylims);
            if thisP<0.001
                startTxt = '***';
            elseif thisP<0.01
                startTxt = '**';
            elseif thisP<0.05
                startTxt = '*';
            end
            text(xs(hi),starY,startTxt,'HorizontalAlignment','center');        
        end
    end
end

set(gca,'XTick',1:nAreas,'XTickLabels',areaLabs);

xlim(xlims);
ylim(ylims);

ylabel(yLab);

legend(hemHs,hemLabs,'Location','SouthEast');


%save figure;
set(gcf,'units','centimeters','pos',[8 8 figWidth figHeight],'color','w');
figTitle = fullfile(p.figures,'LocalizerWordResps_LateralityIndices.eps');
exportfig(gcf,figTitle,'Format','eps','bounds','loose','color','rgb','LockAxes',0,'FontMode','fixed','FontSize',fontSize);


%% Figure S3 for White, Palmer, Boynton & Yeatman, PNAS 2019
% Makes a bar plot of mean adjusted r^2 values from the one- and
% two-channel spatial encoding models, in VWFA-1 and VWFA-2 of both hemispheres. 
%
% by Alex L. White at the University of Washington, 2019
%
function FigS3_EncodingModelAdjustedR2s()

%% set paths
%add whole analysis code directory to the path
analysisDir = fileparts(fileparts(which(mfilename)));
addpath(genpath(analysisDir));

p = getPaths(); 

%% load data
resFile = fullfile(p.results,'AllSubjChannelResponses.mat');
load(resFile);

%% open file to print stats
statsFile = fullfile(p.stats,'StatsS3_ChannelModelFitStats.txt');

statsF = fopen(statsFile,'w');

fprintf(statsF,'ANALYSIS OF VWFA_ATTN5: \nInverted encoding model fit statistics\n');

%% pull out data

valsByIndex = allR.rSqrValsByIndex;

areasToPlot = [find(strcmp(valsByIndex.brainArea,'VWFA_1')) find(strcmp(valsByIndex.brainArea,'VWFA_2'))];
nAreas = length(areasToPlot);
brainAreas = valsByIndex.brainArea(areasToPlot);
areaLabs = cell(1,nAreas);
%take out hyphens
for ai=1:nAreas
    al = brainAreas{ai};
    al(al=='_') = '-';
    areaLabs{ai} = al;
end

hemIs = 1:2;
nHems = length(hemIs);
hemLabs = valsByIndex.hemisphere(hemIs);

modelTypeIs = 1:2; %plot the one-channel model and the 2-channel model
nModelTypes = length(modelTypeIs);
modelLabels = cell(1,nModelTypes);
for mi = 1:nModelTypes
    modelLabels{mi} = sprintf('%i channel', valsByIndex.modelNChannels(mi));
end

rSqrType = 'adjusted'; %plot R2 adjusted for number of parameters in each model type
rSqrTypeI = find(strcmp(valsByIndex.rSqrType, rSqrType));
plotAdjustedR2 = strcmp(rSqrType,'adjusted');

rs = squeeze(allR.rSqrs(areasToPlot,hemIs,modelTypeIs,rSqrTypeI,:));

%% compute means and pairwise differences between model types, and print stats 

fprintf(statsF,'\n\n========================================================\n');
if plotAdjustedR2
    fprintf(statsF,'ADJUSTING R2 FOR NUMBER OF PARAMETERS\n\n');
    ylab = 'adjusted R^2';
    figTitle = fullfile(p.figures,'FigS3_AdjRSqrs_VWFAs.eps');
else
    fprintf(statsF,'**NOT** ADJUSTING R2 FOR NUMBER OF PARAMETERS\n\n');
    ylab = 'R^2';
    figTitle = fullfile(p.figures,'FigS3_RSqrs_VWFAs.eps');
end

mbs = NaN(nAreas, nHems, nModelTypes);
sembs = mbs;
bootSigs = NaN(nAreas, nHems);
bootPs = bootSigs;

for hi=1:nHems
    for ai=1:nAreas
        areaRs = squeeze(rs(ai,hi,:,:));
        %Select subjects with non-nan for both model types
        goodS = all(~isnan(areaRs),1);
        
        areaRs = areaRs(:,goodS);
        mbs(ai,hi,:) = mean(areaRs,ndims(areaRs));
        sembs(ai,hi,:) = standardError(areaRs,ndims(areaRs));
        
        fprintf(statsF,'\n-----------------------------');
        fprintf(statsF,'\n-----------------------------');
        fprintf(statsF,'\nANALYSIS OF %s %s\n', hemLabs{hi},areaLabs{ai});
        
        if plotAdjustedR2
            fprintf(statsF,'\n Mean (SEM) adjusted R2s for two model types:\n');
        else
            fprintf(statsF,'\n Mean (SEM) R2s (NOT adjusted) for two model types:\n');
        end
        
        for modelTypeI = 1:nModelTypes
            modelType = modelTypeIs(modelTypeI);
            fprintf(statsF,'%s:\t %.3f (%.3f)\n', modelLabels{modelTypeI}, mbs(ai,hi,modelType), sembs(ai,hi,modelType));
        end
        
        typeDiffs = diff(areaRs,1,1);
        [diffCI,~,bootDist] = boyntonBootstrap(@nanmean,typeDiffs,5000,95);
        bootSigs(ai,hi) = all(diffCI>0) || all(diffCI<0);
        bootPs(ai,hi) = getBootPs(bootDist',0);
        
        fprintf(statsF,'\n95%% CI of difference: [%.3f %.3f], p=%.4f\n', diffCI(1), diffCI(2),bootPs(ai,hi));
    end
end


%% set up plot 

% bar colors
fillColors = hsv2rgb([0.67 0.4 0.9; 0.67 0.8 0.5]);
edgeColors = fillColors;

semLineWid = 1;
fontSize = 10;
axLineWidth = 1;
bwid = 0.15;

ylims = [0 1.0];
xlims = [0.5 nAreas+0.5];
xjitt = [-0.15 0.15];

dYTick = 0.2;

%% plot

figure;
for hi=1:nHems
    subplot(1,nHems,hi);
    hold on;
    
    for ai = 1:nAreas
        hands = zeros(1,2);
        for modelTypeI = 1:length(modelTypeIs)
            modelType = modelTypeIs(modelTypeI);
            xval = ai + xjitt(modelType);
            bx = xval + [-bwid/2 bwid/2];
            by = [0 mbs(ai,hi,modelType)];
            
            vertx=[bx; bx];
            verty=[by fliplr(by)];
            
            hands(modelTypeI)=fill(vertx(:), verty(:),fillColors(modelType,:),'EdgeColor',edgeColors(modelType,:),'LineWidth',2);
            %do it again in reverse order to avoid diagonal white line
            fill(flipud(vertx(:)), flipud(verty(:)),fillColors(modelType,:),'EdgeColor',edgeColors(modelType,:),'LineWidth',2);
            
            %plot error bar
            plot([xval xval],mbs(ai,hi,modelType)+[-1 1]*sembs(ai,hi,modelType),'k-','LineWidth',semLineWid);
        end
        
        %significance of comparing model types
        if bootSigs(ai,hi)
            textx =  ai; % + xjitt(modelTypesToCompare(compI,1));
            text(textx,0.75,'*','Color',[0 0 0],'HorizontalAlignment','center');
        end
    end
    
    xlim(xlims);
    ylim(ylims);
    set(gca,'XTick',1:nAreas,'YTick',ylims(1):dYTick:ylims(2),'LineWidth',axLineWidth);
    set(gca,'XTickLabels',areaLabs);
    if hi==1
        ylabel(ylab)
        legend(hands,modelLabels);
    else
        set(gca,'YTickLabels',{});
    end
    
    title([hemLabs{hi} ' Hemisphere']);  
end

%save figure;
set(gcf,'units','centimeters','pos',[8 8 11.4 6],'color','w');
exportfig(gcf,figTitle,'Format','eps','bounds','loose','color','rgb','LockAxes',0,'FontMode','fixed','FontSize',fontSize);

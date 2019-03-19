%% Stats to accompany Figures 3 and S4 in White, Palmer, Boynton & Yeatman, PNAS 2019
% Analysis of channel responses in both hemisphere VWFAs.
% by Alex L. White at the University of Washington, 2019
%
function Stats3_ChannelRespLME()

%% set paths
%add whole analysis code directory to the path
analysisDir = fileparts(fileparts(which(mfilename)));
addpath(genpath(analysisDir));

p = getPaths(); 

%% load data
resFile = fullfile(p.results,'AllSubjChannelResponses.mat');
load(resFile);

%% Open a file to print stats to:
statsFile = fullfile(p.stats,'Stats3_VWFAChannelResponses.txt');
diary(statsFile);
statsF = fopen(statsFile,'w');

fprintf(1,'Analys of channel responses from inverted encoding model in VWFAs\n');
fprintf(1,'Unless otherwise stated, ''DummyVarCoding'' = ''effects'', so we can estimate main effects and interactions:\n');
fprintf(1,'Recoding the cue as: distributed; focal *cued*; or focal *uncued*');

%% Extract channel responses
areasToPlot = [find(strcmp(allR.valsByIndex.brainArea,'VWFA_1')) find(strcmp(allR.valsByIndex.brainArea,'VWFA_2'))];
nAreas = length(areasToPlot);
brainAreas = allR.valsByIndex.brainArea(areasToPlot);
%replace underscores with hyphens in brain area names 
areaLabs = cell(1,nAreas);
for ai=1:nAreas
    al = brainAreas{ai};
    al(al=='_') = '-';
    areaLabs{ai} = al;
end

hemIs = 1:2;
nHems = length(hemIs);
hemLabs = allR.valsByIndex.hemisphere(hemIs);

channIs = 1:2;
channelLabels = allR.valsByIndex.channel(channIs);
nChann = length(channIs);

condsToPlot = 1:3;
condLabs = allR.valsByIndex.condition(condsToPlot);
nCond = numel(condLabs);

%individual subject chanenl responses: 
rs = squeeze(allR.twoChannelResponses(areasToPlot,hemIs,channIs,condsToPlot,:));

nSubj = size(rs,ndims(rs));

%% create a table
%make a table
nRows      = numel(rs);
subject    = NaN(nRows,1);
region     = cell(nRows,1);
hemisphere = cell(nRows,1);
channel    = cell(nRows,1);
cue        = cell(nRows,1);
resp       = NaN(nRows,1);

r = 0;
for si = 1:size(rs,ndims(rs))
    for ai = 1:nAreas
        for hi = 1:nHems
            for ki = 1:nChann
                for ci = 1:nCond
                    
                    r = r+1;
                    subject(r) = si;
                    
                    region(r) = areaLabs(ai);
                    hemisphere(r) = hemLabs(hi);
                    channel(r) = channelLabels(ki);
                                        
                    if strcmp(condLabs(ci),'Distributed cue')
                        cue(r) = {'distributedCue'};
                        
                    elseif strcmp(condLabs(ci), 'Focal cue left')
                        if strcmp(channelLabels(ki),'Left')
                            %left channel is "cued"
                            cue(r) = {'focalCued'};
                        elseif strcmp(channelLabels(ki),'Right')
                            %rigth channel is "uncued"
                            cue(r) = {'focalUncued'};
                        end
                        
                    elseif strcmp(condLabs(ci), 'Focal cue right')
                        if strcmp(channelLabels(ki),'Left')
                            %left channel is "uncued"
                            cue(r) = {'focalUncued'};
                        elseif strcmp(channelLabels(ki),'Right')
                            %right channel is "cued" 
                            cue(r) = {'focalCued'};
                        end  
                    end
                    
                    resp(r) = rs(ai,hi,ki,ci,si);
                end
            end
        end
    end
end

subject = categorical(subject);
d = table(subject, region, hemisphere, channel, cue, resp);

%% Run LMEs: effect of channel (left vs right), cue, and the interaction of channel of cue, in each ROI
for hem = 1:2
    hemSub = strcmp(d.hemisphere, hemLabs{hem});
    
    %SEPARATE ANALYSIS OF EACH ROI
    for ai=1:numel(areaLabs)
        thisArea = areaLabs{ai};
        areaSub = strcmp(d.region,thisArea);
        fprintf(1,'\n------------------------------------------------------\n');
        fprintf(1,'------------------------------------------------------\n');
        fprintf(1,'ANAYSLIS OF %s %s ONLY', hemLabs{hem},thisArea);
        fprintf(1,'\n------------------------------------------------------\n');
        fprintf(1,'------------------------------------------------------\n');
        
        %find subjects who have at least some data from this region
        goodSubj = false(size(d,1),1);
        subjsToInclude = [];
        for sii=1:nSubj
            subjResps = d.resp(hemSub & areaSub & d.subject==categorical(sii));
            isGood = ~all(isnan(subjResps));
            goodSubj(d.subject==categorical(sii)) = isGood;
            if isGood, subjsToInclude = [subjsToInclude sii]; end
        end
                
        %% Selective attention effects: focal cued vs focal uncued responses, in each channel 
        
        %1. Two-way effects of channel and cue, excluding distributed cue, with random slopes for each subject
        fprintf(1,'\nLME for two-way effects of channel (2), and cue (2), excluding distributed-cue, *with random slopes for each subject*, in %s %s\n', hemLabs{hem},thisArea);

        subst = hemSub & areaSub & goodSubj & ~strcmp(d.cue,'distributedCue');

        eqtn = 'resp ~ channel*cue + (channel*cue | subject)';
        lme1 = fitlme(d(subst,:), eqtn, 'DummyVarCoding','effects');
        lme1.disp;
        fprintf(1,'\n---Corresponding ANOVA---\n');
        display(lme1.anova);

        %1b: Test the simpler model to see if including random slopes was worth justified: 
        eqtn = 'resp ~ channel*cue + (1 | subject)';
        lme1b = fitlme(d(subst,:), eqtn, 'DummyVarCoding','effects');
        fprintf(1,'\n---Corresponding ANOVA---\n');
        
        fprintf(1,'\nComparing models with and without random slopes by subject\n');
        compare(lme1b,lme1)
        
        %Compute mean and SEM selective effect sizes 
        subTable  = d(subst,:); 
        chanNames = unique(subTable.channel);
        cueNames = unique(subTable.cue);
        crs = NaN(2,2,nSubj);
        for chni = 1:2
            for cui=1:2
                condSub = strcmp(subTable.channel,chanNames{chni}) & strcmp(subTable.cue,cueNames{cui});
                crs(chni,cui,subTable.subject(condSub)) = subTable.resp(condSub,:);
            end
        end
        
        selectEs = -1*squeeze(diff(crs,1,2));
        meanSelectEs = nanmean(selectEs,2); 
        semSelectEs = standardError(selectEs,2);
        
        selectEsChnMean = squeeze(nanmean(selectEs,1));
        meanSelectEChnMean = nanmean(selectEsChnMean); 
        semSelectEChnMean = standardError(selectEsChnMean);
        
        fprintf(1,'\nMean Selective Attention Effects (%s - %s) in %s %s:\n', cueNames{1}, cueNames{2}, hemLabs{hem},thisArea);
        for chni = 1:2
            fprintf(1,'%s channel: %.3f +/- %.3f\n', chanNames{chni}, meanSelectEs(chni), semSelectEs(chni));
        end
        fprintf(1,'Avg across channels: %.3f +/- %.3f\n', meanSelectEChnMean, semSelectEChnMean);

        %Compute ratio cued/uncued of mean channel responses
        meanCRs = nanmean(crs,3); 
        meanSelectRs = meanCRs(:,1)./meanCRs(:,2);
        meanSelectRsChnMean = nanmean(meanSelectRs);
        
        fprintf(1,'\nMean Selective Attention RATIOS (%s/%s) in %s %s:\n*computed by averaging channel responses over subjects first*\n', cueNames{1},cueNames{2},hemLabs{hem},thisArea);
        for chni = 1:2
            fprintf(1,'%s channel: %.3f\n', chanNames{chni}, meanSelectRs(chni));
        end
        fprintf(1,'Avg across channels: %.3f\n', meanSelectRsChnMean);

        %% Divided attention effects: focal cued vs distributed cue 
        
        %2. Two-way effects of channel and cue, excluding focal uncued
        %(divided attention effects), with random slopes by subject
        subst = hemSub & areaSub & goodSubj & ~strcmp(d.cue,'focalUncued');
        fprintf(1,'\nLME for two-way effects of channel (2), and cue (2), excluding focal uncued, *with random slopes for each subject*, in %s %s\n', hemLabs{hem},thisArea);
        eqtn = 'resp ~ channel*cue + (channel*cue | subject)';
        lme2 = fitlme(d(subst,:), eqtn, 'DummyVarCoding','effects');
        lme2.disp;
        fprintf(1,'\n---Corresponding ANOVA---\n');
        display(lme2.anova);
        
        fprintf(1,'\nComparing the last two models with and without random effects by subject\n');
        eqtn = 'resp ~ channel*cue + (1 | subject)';
        lme2b = fitlme(d(subst,:), eqtn, 'DummyVarCoding','effects');   
        compare(lme2b,lme2)
                  
        %Compute mean and SEM divided effect sizes 
        subTable  = d(subst,:); 
        
        chanNames = unique(subTable.channel);
        cueNames = unique(subTable.cue);
        crs = NaN(2,2,nSubj);
        for chni = 1:2
            for cui=1:2
                condSub = strcmp(subTable.channel,chanNames{chni}) & strcmp(subTable.cue,cueNames{cui});
                crs(chni,cui,subTable.subject(condSub)) = subTable.resp(condSub,:);
            end
        end
        
        dividEs = -1*squeeze(diff(crs,1,2));
        meanDividEs = nanmean(dividEs,2); 
        semDividEs = standardError(dividEs,2);
        
        dividEsChnMean = squeeze(nanmean(dividEs,1));
        meanDividEChnMean = nanmean(dividEsChnMean); 
        semDividEChnMean = standardError(dividEsChnMean);
        
        fprintf(1,'\nMean Divided Attention Effects (%s-%s) in %s %s:\n', cueNames{1}, cueNames{2}, hemLabs{hem},thisArea);
        for chni = 1:2
            fprintf(1,'%s channel: %.3f +/- %.3f\n', chanNames{chni}, meanDividEs(chni), semDividEs(chni));
        end
        fprintf(1,'Avg across channels: %.3f +/- %.3f\n', meanDividEChnMean, semDividEChnMean);

        %Compute ratio distrubted/focal of mean channel responses
        meanCRs = nanmean(crs,3); 
        meanDividRs = meanCRs(:,1)./meanCRs(:,2);
        meanDividRsChnMean = nanmean(meanDividRs);
        
        fprintf(1,'\nMean Divided Attention RATIOS (%s/%s) in %s %s:\n*computed by averaging channel repsonses over subjects first*\n', cueNames{1}, cueNames{2}, hemLabs{hem},thisArea);
        for chni = 1:2
            fprintf(1,'%s channel: %.3f\n', chanNames{chni}, meanDividRs(chni));
        end
        fprintf(1,'Avg across channels: %.3f\n', meanDividRsChnMean);
        
        %% finally, compute mean difference across channels 
        subst = hemSub & areaSub & goodSubj;
        subTable  = d(subst,:); 
        
        chanNames = unique(subTable.channel);
        cueNames = unique(subTable.cue);
        crs = NaN(2,3,nSubj);
        for chni = 1:2
            for cui=1:3
                condSub = strcmp(subTable.channel,chanNames{chni}) & strcmp(subTable.cue,cueNames{cui});
                crs(chni,cui,subTable.subject(condSub)) = subTable.resp(condSub,:);
            end
        end
        condMeans = squeeze(nanmean(crs,2));
        channDiffs = squeeze(diff(condMeans,1,1));
        meanChannelDiff = nanmean(channDiffs); 
        semChannelDiff = standardError(channDiffs);
        fprintf(1,'\n\nAveraging over cue conds, mean (sem) difference in %s-%s channels: %.3f (%.3f)\n', chanNames{2}, chanNames{1}, meanChannelDiff, semChannelDiff);
    end
end


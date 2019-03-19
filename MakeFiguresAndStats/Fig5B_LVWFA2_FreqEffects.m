%% Figure 5B for White, Palmer, Boynton & Yeatman, PNAS 2019
% Makes a bar plot of the effects of each word's lexical frequency, in each
% cue condition, on mean BOLD response in left VWFA-2. 
%
% by Alex L. White at the University of Washington, 2019
function Fig5B_LVWFA2_FreqEffects(figH)

%% set paths 
p = getPaths(); 

%% load data
resF = fullfile(p.results, 'AllSubjWordFreqResponses.mat');
load(resF); 

%% open a stats file 
statsFile = fullfile(p.stats,'Stats5B_FreqEffectsMeanBOLD.txt'); 
statsF = fopen(statsFile,'w'); 
diary(statsFile); 

fprintf(1,'VWFA-Attn5: Results of trial-wise event-related analysis, coding for the lexical frequency of each word\n');

%% pull out data

areasToAnalyze = [find(strcmp(allR.valsByIndex.brainArea,'VWFA_1')) find(strcmp(allR.valsByIndex.brainArea,'VWFA_2'))];
nAreas = length(areasToAnalyze);
brainAreas = allR.valsByIndex.brainArea(areasToAnalyze);

areaLabs = cell(1,nAreas);
%take out hyphens
for ai=1:nAreas
    al = brainAreas{ai};
    al(al=='_') = '-';
    areaLabs{ai} = al;
end

hemIs = 1:2;
nHems = length(hemIs);
hemLabs = allR.valsByIndex.hemisphere(hemIs);

% parse each condition into three variables: cue condition, left word
% frequency, right word frequence
% and give each condition 3 numbers to identify it 
cns = allR.valsByIndex.condition;
nReg = length(cns);
condNames = cell(nReg,3);
%parse each condition name by the underscores 
for cni=1:nReg
    clear uis;
    cn = cns{cni};
    uis = find(cn=='_');
    uis = [0 uis length(cn)+1];
    for varI = 1:(length(uis)-1)
        condNames{cni,varI} = cn((1+uis(varI)):(uis(varI+1)-1));
    end
end

%rename
condNames(strcmp(condNames,'distrib')) = {'Distributed'};
condNames(strcmp(condNames,'focalLeft')) = {'Focal Left'};
condNames(strcmp(condNames,'focalRight')) = {'Focal Right'};
condNames(strcmp(condNames,'leftLow')) = {'L Low'};
condNames(strcmp(condNames,'leftHigh')) = {'L High'};
condNames(strcmp(condNames,'rightLow')) = {'R Low'};
condNames(strcmp(condNames,'rightHigh')) = {'R High'};

uniqueCondNames = cell(1,3);
nCondsByVar = NaN(1,3);
variableNames = cell(1,3);
for varI = 1:3
    uniqueCondNames(varI) = {unique(condNames(:,varI))};
    nCondsByVar(varI) = numel(uniqueCondNames{varI});
    if any(strcmp(uniqueCondNames{varI},'Distributed'))
        variableNames{varI} = 'Cue';
    else
        firstLetts = cell(1,nCondsByVar(varI));
        for fli=1:nCondsByVar(varI)
            firstLetts{fli} = uniqueCondNames{varI}{fli}(1);
        end
        if all(strcmp(firstLetts,'L'))
            variableNames{varI} = 'Left freq';
        elseif all(strcmp(firstLetts,'R'))
            variableNames{varI} = 'Right freq';
        end
    end
end
%give each condition a number
condNums = NaN(nReg,3);
for cni=1:length(cns)
    for varI=1:3
        condNums(cni,varI) = find(strcmp(condNames{cni,varI},uniqueCondNames{varI}));
    end
end

%% Run statistics on each area

nSubj = size(allR.meanBetas, ndims(allR.meanBetas));

% also reorganize all the data into a matrix: area x hemisphere x channel x
% cue condition x left feature x right feature
datMat = NaN(nAreas, nHems, nCondsByVar(1), nCondsByVar(2), nCondsByVar(3), nSubj);

for aii=1:nAreas
    ai = areasToAnalyze(aii);
    for hii=1:nHems
        hi = hemIs(hii);
        for v1i = 1:nCondsByVar(1)
            for v2i = 1:nCondsByVar(2)
                for v3i = 1:nCondsByVar(3)
                    condI = condNums(:,1)==v1i & condNums(:,2)==v2i & condNums(:,3)==v3i;
                    theseDats = squeeze(allR.meanBetas(ai,hi,condI,:));
                    datMat(aii,hii,v1i,v2i,v3i,:) = theseDats;
                end
            end
        end
        
    end
end


%compute effects of variable2 (left freq), averaging
%over variable3 (right freq)
v3Avg = squeeze(nanmean(datMat,5));
v2Effects = squeeze(v3Avg(:,:,:,1,:) - v3Avg(:,:,:,nCondsByVar(2),:)); %first minus last
meanEffectsOfVar2 = nanmean(v2Effects,ndims(v2Effects));
semEffectsOfVar2 = standardError(v2Effects,ndims(v2Effects));

%compute effects of variable3 (right freq), averaging
%over variable2 (left freq)
v2Avg = squeeze(nanmean(datMat,4));
v3Effects = squeeze(v2Avg(:,:,:,1,:) - v2Avg(:,:,:,nCondsByVar(3),:)); %first minus last
meanEffectsOfVar3 = nanmean(v3Effects,ndims(v3Effects));
semEffectsOfVar3 = standardError(v3Effects,ndims(v3Effects));

var2effectName = sprintf('%s - %s', uniqueCondNames{2}{1}, uniqueCondNames{2}{end});
var3effectName = sprintf('%s - %s', uniqueCondNames{3}{1}, uniqueCondNames{3}{end});


%Run stats on effects!  
effectsOfVar2_95CIs = NaN([size(meanEffectsOfVar2) 2]);
effectsOfVar3_95CIs = NaN([size(meanEffectsOfVar3) 2]);
var2EffectBootPs = NaN(size(meanEffectsOfVar2));
var3EffectBootPs = NaN(size(meanEffectsOfVar3));


for aii=1:size(v2Effects,1)
    for hii=1:size(v2Effects,2)
        
        %make a table
        cues = {};
        effects = [];
        subjs = [];
        sides = {};

        fprintf(1,'\n\n===============================================\n');
        fprintf(1,'ANALYSIS OF EFFECTS OF FREQUENCY  IN %s %s:\n\n', hemLabs{hii}, areaLabs{aii});
     
        
        for v1i=1:size(v2Effects,3)
            
            
            theseE2s = squeeze(v2Effects(aii,hii,v1i,:));
            [e2CI, ~, e2Dist] = boyntonBootstrap(@nanmean,theseE2s,1000,95,true);
            effectsOfVar2_95CIs(aii,hii,v1i,:) = e2CI;
            var2EffectBootPs(aii,hii,v1i) = getBootPs(e2Dist',0);
            
            theseE3s = squeeze(v3Effects(aii,hii,v1i,:));
            [e3CI, ~, e3Dist] = boyntonBootstrap(@nanmean,theseE3s,1000,95,true);
            effectsOfVar3_95CIs(aii,hii,v1i,:) = e3CI;
            var3EffectBootPs(aii,hii,v1i) = getBootPs(e3Dist',0);
            
            %print stats 
            fprintf(1,'\nEffect sizes in cue condition %s:\n', uniqueCondNames{1}{v1i});
            fprintf(1,'%s: mean=%.4f, SEM = %.3f, 95%%CI = [%.3f %.3f], bootP=%.5f,', var2effectName, nanmean(theseE2s), standardError(theseE2s'), e2CI(1), e2CI(2), var2EffectBootPs(aii,hii,v1i))
            [~,tp,~,tst] = ttest(theseE2s);
            fprintf(1,'\n\tt(%i)=%.3f,p=%.5f\n', tst.df, tst.tstat, tp);
            
            fprintf(1,'%s: mean=%.4f, SEM = %.3f, 95%%CI = [%.3f %.3f], bootP=%.5f,', var3effectName, nanmean(theseE3s), standardError(theseE3s'), e3CI(1), e3CI(2), var3EffectBootPs(aii,hii,v1i))
            [~,tp,~,tst] = ttest(theseE3s);
            fprintf(1,'\n\tt(%i)=%.3f,p=%.5f\n', tst.df, tst.tstat, tp);
            
            
            %assemble data for table
            cues = cat(1,cues,uniqueCondNames{1}(ones(length(theseE2s)+length(theseE3s),1)*v1i));
            effects = [effects; theseE2s; theseE3s]; 
            sideLabels = cat(1,repmat({'left'},size(theseE2s)),repmat({'right'},size(theseE3s)));
            sides = cat(1,sides,sideLabels);

            subjs = [subjs; (1:nSubj)'; (1:nSubj)'];

        end
        
        T2 = table; 
        T2.subject = subjs;
        T2.effect = effects; 
        T2.side = sides; 
        T2.cue = cues;
        
        fprintf(1,'\n2way analysis of cue (all conds) and word side on the effect of frequency (%s %s)', hemLabs{hii}, areaLabs{aii});
        
        eqtn = 'effect ~ cue*side'; 
        lme = fitlme(T2,eqtn,'DummyVarCoding','effect');
        display(lme);
        display(lme.anova);
        
        %doing this with random effects by subject leads to the same result 
        
        fprintf(1,'\n2way analysis of cue (JUST FOCAL!) and word side on the effect of frequency (%s %s)', hemLabs{hii}, areaLabs{aii});
        
        cueSubst = ~strcmp(T2.cue,'Distributed');
        
        lme = fitlme(T2(cueSubst,:),eqtn,'DummyVarCoding','effect');
        display(lme);
        display(lme.anova);
        
    end
end

diary off;


%% set up plot

cueHues = [0.8 0.6 0.4];
cueSats = [0.9 0.9 0.9];
cueVals = [0.7 0.7 0.7];
cueColrs = hsv2rgb([cueHues' cueSats' cueVals']);


edgeColors = NaN(size(cueColrs,1),2,3);
edgeColors(:,1,:) = cueColrs;
edgeColors(:,2,:) = cueColrs;

fillColors = edgeColors;
fillColors(:,1,:) = 1;

opt.edgeColors = edgeColors;
opt.fillColors = fillColors;
opt.errorBarColors = zeros(size(opt.edgeColors));

opt.barWidth = 0.09;
opt.edgeLineWidth = 1.5;
opt.errorBarWidth = 1;
opt.level1Sep = 0.25;
opt.level2Sep = 0.15;
opt.xAxisMargin = 0.25;

opt.ylims = [-0.04 0.04];
opt.yLab = '\Delta psc';

opt.xTickLabs = {}; %uniqueCondNames{1};

opt.legendLabs = {var2effectName,var3effectName};
opt.legendLoc = 'NorthWest';
opt.lev1ForLegend = 2;


errorBarType = 1; %1=SEM, 2=CI

%% plot
figure(figH);

%plot frequency effects of each word, in each cue condition, in just one area, one hemisphere: 
areaNameToPlot = 'VWFA_2';
areasToPlot = find(strcmp(brainAreas,areaNameToPlot));
nAreas = length(areasToPlot); 
hemsToPlot = find(strcmp(hemLabs,'Left'));
nHems = length(hemsToPlot);

nRows = 1;
nCols = 2;
subi = 2; %fixed to larger figure
for aii=1:nAreas
    ai = areasToPlot(aii);
    for hii=1:nHems
        hi = hemsToPlot(hii);
        
        subplot(nRows,nCols,subi); hold on;
        
        m2 = squeeze(meanEffectsOfVar2(ai,hi,:));
        m3 = squeeze(meanEffectsOfVar3(ai,hi,:));
        ms = [m2 m3];
        
        if errorBarType==1
            %SEMs:
            e2 = squeeze(semEffectsOfVar2(ai,hi,:));
            e3 = squeeze(semEffectsOfVar3(ai,hi,:));
            es = [e2 e3];
        else
            %95% CIs
            e2 = squeeze(effectsOfVar2_95CIs(ai,hi,:,:));
            e3 = squeeze(effectsOfVar3_95CIs(ai,hi,:,:));
            
            es = NaN([size(ms) 2]);
            for ebi=1:2
                es(:,1,ebi) = e2(:,ebi);
                es(:,2,ebi) = e3(:,ebi);
            end
        end
        
        opt.doLegend = false; %aii==1 & hii==1;
        
        barCenters = barPlot_AW(ms, es, opt);
        title([hemLabs{hii} ' ' areaLabs{ai}]);
        
        set(gca,'TitleFontSizeMultiplier',1,'TitleFontWeight','normal');
        
        if hii>1
            set(gca,'YTickLabel',{});
            ylabel('');
        end
      
        %add stars for significance
        for v1i = 1:size(ms,1)
            ps = [var2EffectBootPs(ai,hi,v1i) var3EffectBootPs(ai,hi,v1i)];
            for jj=1:2
                if ps(jj)<0.05
                    starY = opt.ylims(1)+0.07*diff(opt.ylims);
                    if ps(jj)<0.001
                        startTxt = '***';
                    elseif ps(jj)<0.01
                        startTxt = '**';
                    elseif ps(jj)<0.05
                        startTxt = '*';
                    end
                    text(barCenters(v1i,jj),starY,startTxt,'HorizontalAlignment','center');
                end
            end
        end
    end
end


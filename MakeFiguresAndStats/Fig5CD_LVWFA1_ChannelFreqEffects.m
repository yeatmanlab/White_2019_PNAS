function Fig5CD_LVWFA1_ChannelFreqEffects(figH)

%% set paths
p = getPaths(); 

%% load data
resF = fullfile(p.data, 'AllSubjChannelResponses_WordFreq.mat');
load(resF);


%% open a stats file
statsFile = fullfile(p.stats,'Stats5C_FreqEffectsChannelResps.txt');
statsF = fopen(statsFile,'w');
diary(statsFile);

fprintf(1,'Effects of lexical frequency on channel responses in VWFAs\n');

%% pull out data

areaToPlot = find(strcmp(allR.valsByIndex.brainArea,'VWFA_1'));
brainArea = allR.valsByIndex.brainArea{areaToPlot};
areaLab = brainArea;
areaLab(areaLab=='_') = '-';

hemI = 1;
hemLab = allR.valsByIndex.hemisphere{hemI};

channIs = 1:2;
nChann = length(channIs);
channLabs = allR.valsByIndex.channel;
channLabs(strcmp(channLabs,'Left')) = {'L. Channel'};
channLabs(strcmp(channLabs,'Right')) = {'R. Channel'};

% parse each condition into three variables: cue condition, left word
% frequency, right word frequence
% and give each condition 3 numbers to identify it 
cns = allR.valsByIndex.condition';
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

nSubj = size(allR.twoChannelResponses, ndims(allR.twoChannelResponses));

%% Run stats 

% reorganize all the data, for just one area into a matrix:  channel x cue
% condition x left word freq x right word freq
datMat = NaN(nChann, nCondsByVar(1), nCondsByVar(2), nCondsByVar(3), nSubj);
for kii=1:nChann
    ki = channIs(kii);
    for v1i = 1:nCondsByVar(1)
        for v2i = 1:nCondsByVar(2)
            for v3i = 1:nCondsByVar(3)
                condI = condNums(:,1)==v1i & condNums(:,2)==v2i & condNums(:,3)==v3i;
                theseDats = squeeze(allR.twoChannelResponses(areaToPlot,hemI,ki,condI,:));
                datMat(kii,v1i,v2i,v3i,:) = theseDats;
                  
            end
        end
    end
end

%compute effects of variable2 (left freq), averaging
%over variable3 (right freq)
v3Avg = squeeze(nanmean(datMat,4));
v2Effects = squeeze(v3Avg(:,:,1,:) - v3Avg(:,:,nCondsByVar(2),:)); %first minus last
meanEffectsOfVar2 = nanmean(v2Effects,ndims(v2Effects));
semEffectsOfVar2 = standardError(v2Effects,ndims(v2Effects));

%compute effects of variable3 (right freq), averaging
%over variable2 (left freq)
v2Avg = squeeze(nanmean(datMat,3));
v3Effects = squeeze(v2Avg(:,:,1,:) - v2Avg(:,:,nCondsByVar(3),:)); %first minus last
meanEffectsOfVar3 = nanmean(v3Effects,ndims(v3Effects));
semEffectsOfVar3 = standardError(v3Effects,ndims(v3Effects));

var2effectName = sprintf('%s-%s', uniqueCondNames{2}{1}, uniqueCondNames{2}{end});
var3effectName = sprintf('%s-%s', uniqueCondNames{3}{1}, uniqueCondNames{3}{end});

%Run stats on effects!
effectsOfVar2_95CIs = NaN([size(meanEffectsOfVar2) 2]);
effectsOfVar3_95CIs = NaN([size(meanEffectsOfVar3) 2]);
var2EffectBootPs = NaN(size(meanEffectsOfVar2));
var3EffectBootPs = NaN(size(meanEffectsOfVar3));


for kii=1:size(v2Effects,1)
    %make a table
    cues = {};
    effects = [];
    subjs = [];
    sides = {};
    
    fprintf(1,'\n\n===============================================\n');
    fprintf(1,'ANALYSIS OF EFFECTS OF frequency IN %s of %s %s:\n\n', channLabs{kii},  hemLab, areaLab);
    
    
    for v1i=1:size(v2Effects,2)
        
        
        theseE2s = squeeze(v2Effects(kii,v1i,:));
        [e2CI, ~, e2Dist] = boyntonBootstrap(@nanmean,theseE2s,1000,95,true);
        effectsOfVar2_95CIs(kii,v1i,:) = e2CI;
        var2EffectBootPs(kii,v1i) = getBootPs(e2Dist',0);
        
        theseE3s = squeeze(v3Effects(kii,v1i,:));
        [e3CI, ~, e3Dist] = boyntonBootstrap(@nanmean,theseE3s,1000,95,true);
        effectsOfVar3_95CIs(kii,v1i,:) = e3CI;
        var3EffectBootPs(kii,v1i) = getBootPs(e3Dist',0);
        
        %print stats
        fprintf(1,'\nEffect sizes in cue condition %s:\n', uniqueCondNames{1}{v1i});
        fprintf(1,'%s: mean=%.4f, SEM = %.3f, 95%%CI = [%.3f %.3f], bootP=%.5f,', var2effectName, nanmean(theseE2s), standardError(theseE2s'), e2CI(1), e2CI(2), var2EffectBootPs(kii,v1i))
        [~,tp,~,tst] = ttest(theseE2s);
        fprintf(1,'\n\tt(%i)=%.3f,p=%.5f\n', tst.df, tst.tstat, tp);
        
        fprintf(1,'%s: mean=%.4f, SEM = %.3f, 95%%CI = [%.3f %.3f], bootP=%.5f,', var3effectName, nanmean(theseE3s), standardError(theseE3s'), e3CI(1), e3CI(2), var3EffectBootPs(kii,v1i))
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
    
    fprintf(1,'\n2way analysis of cue (all conds) and word side on the effect of frequency (%s, %s %s)', channLabs{kii}, hemLab, areaLab);
    
    eqtn = 'effect ~ cue*side';
    lme = fitlme(T2,eqtn,'DummyVarCoding','effect');
    display(lme);
    display(lme.anova);
    
    %doing this with random effects by subject leads to the same result
    
    fprintf(1,'\n2way analysis of cue (JUST FOCAL!) and word side on the effect of frequency (%s, %s %s)',  channLabs{kii}, hemLab, areaLab);
    
    cueSubst = ~strcmp(T2.cue,'Distributed');
    
    lme = fitlme(T2(cueSubst,:),eqtn,'DummyVarCoding','effect');
    display(lme);
    display(lme.anova);
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

opt.ylims = [-0.3 0.3];
opt.yLab = '\Delta response';


opt.xTickLabs = uniqueCondNames{1};
opt.xTickLabs{strcmp(opt.xTickLabs,'Distributed')} = 'Distr.';
opt.xTickLabs{strcmp(opt.xTickLabs,'Focal Left')} = 'Foc.L';
opt.xTickLabs{strcmp(opt.xTickLabs,'Focal Right')} = 'Foc.R';

opt.legendLabs = {var2effectName,var3effectName};
opt.legendLoc = 'NorthWest';
opt.lev1ForLegend = 2;

%% plot effects

nRows = 2;
nCols = 2;

errorBarType = 1; %1=SEM, 2=CI

figure(figH);
subi=2; %fixed to larger figure

for kii=1:nChann
    
    subi=subi+1;
    subplot(nRows,nCols,subi); hold on;
    
    m2 = squeeze(meanEffectsOfVar2(kii,:))';
    m3 = squeeze(meanEffectsOfVar3(kii,:))';
    ms = [m2 m3];
    
    if errorBarType==1
        %SEMs:
        e2 = squeeze(semEffectsOfVar2(kii,:))';
        e3 = squeeze(semEffectsOfVar3(kii,:))';
        es = [e2 e3];
    else
        %95% CIs
        e2 = squeeze(effectsOfVar2_95CIs(kii,:,:));
        e3 = squeeze(effectsOfVar3_95CIs(kii,:,:));
        
        es = NaN([size(ms) 2]);
        for ebi=1:2
            es(:,1,ebi) = e2(:,ebi);
            es(:,2,ebi) = e3(:,ebi);
        end
    end
    
    opt.doLegend = false;  
    
    barCenters = barPlot_AW(ms, es, opt);
    
    title([hemLab ' ' areaLab ' ' channLabs{kii}]);
    set(gca,'TitleFontSizeMultiplier',1,'TitleFontWeight','normal');

    if kii>1
        ylabel('');
    end
    
    %add stars for significance
    for v1i = 1:size(ms,1)
        ps = [var2EffectBootPs(kii,v1i) var3EffectBootPs(kii,v1i)];
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
%% Figure 5A for White, Palmer, Boynton & Yeatman, PNAS 2019
% Makes a bar plot of the effects of each word's lexical frequency, in each
% cue condition, on behavioral accuracy expressed as d'. 
%
% by Alex L. White at the University of Washington, 2019

%% plot effect of target frequency and cue condition for each target side
function Fig5A_DprimeSidexFrequency(figH)

%% set paths
%add whole analysis code directory to the path
analysisDir = fileparts(fileparts(which(mfilename)));
addpath(genpath(analysisDir));

p = getPaths(); 
%% load data
resFile = fullfile(p.data,'AllSubjBehavior_WordFreq.mat');
load(resFile);

allDs = allR.dprime;
allUDs = allR.uncuedDprime;

%% extract data 

%dprime for cued stimulus
conds.cueConds = [2 3];
conds.targSide = 2:3;
conds.targLength = 1;
conds.distLength = 1;
conds.targFreqBin = 2:length(rAvg.valsByIndex.targFreqBin);
conds.distFreqBin = 1;
conds.targPres = 1;

cueLabs = {'Distributed','Focal'};
sideLabs = {'Left','Right'};
freqLabs = {'Low','High'};

ds = squeeze(allDs(conds.cueConds, conds.targSide, conds.targLength, conds.distLength, conds.targFreqBin, conds.distFreqBin, conds.targPres,:));

%dprime for *uncued* stimulus
condsU.cueConds = [2 3];
condsU.targSide = 2:3;
condsU.targLength = 1;
condsU.distLength = 1;
condsU.targFreqBin = 1;
condsU.distFreqBin = 2:length(rAvg.valsByIndex.distFreqBin);
condsU.targPres = 1;

us = squeeze(allUDs(condsU.cueConds, condsU.targSide, condsU.targLength, condsU.distLength, condsU.targFreqBin, condsU.distFreqBin, condsU.targPres,:));

%% run some stats!
statsFile = fullfile(p.stats,'Stats5A_AccuracyByFrequency.txt');
diary(statsFile);

statsF = fopen(statsFile,'w');

fprintf(1,'Analysis of d'' in the scanner as a function of cue condition, word side, and lexical frequency bin\n');

%make a table
nRows = numel(ds);

dprimes = NaN(nRows,1);
cueConds = cell(nRows,1);
freqBins = cell(nRows,1);
sides = cell(nRows,1);
subjects = NaN(nRows,1);
rowi=0;
for ci=1:length(conds.cueConds)
    for tsi=1:length(conds.targSide)
        for tfi=1:length(conds.targFreqBin)
            for si=1:rAvg.nsubj
                rowi=rowi+1;
                dprimes(rowi) = ds(ci,tsi,tfi,si);
                cueConds{rowi} = cueLabs{ci};
                sides{rowi} = sideLabs{tsi};
                freqBins{rowi} = freqLabs{tfi};
                subjects(rowi) = si;
            end
        end
    end
end

T = table;
T.subject = categorical(subjects);
T.side = sides;
T.cueCond = cueConds;
T.freqBin = freqBins;
T.dprime = dprimes;

fprintf(1,'LME on dprime for 3-way effects of cue, side, and frequency bin, with random effects by subject\n');

eqtn = 'dprime ~ cueCond*side*freqBin + (1|subject)';

lme = fitlme(T,eqtn,'DummyVarCoding','effects');
display(lme);
display(lme.anova)

%% compute EFFECTS of each word's frequency in each cue condition
targEffects = squeeze(diff(ds,1,3));  %effect of target word's frequency on d' for target stimulus
distEffects = squeeze(diff(us,1,3));  %effect of distractor word's frequency on d' for distractor stimulus
allEffects = NaN(3,2,rAvg.nsubj); %cue (D, FL, FR) x word side (L,R)

xLabs = cell(1,3);
for cueSide = 1:3
    switch cueSide
        case 1 %distributed cue
            cueRow = 1;
            sideCol = 1:2;
            allEffects(cueSide,:,:) = targEffects(cueRow, sideCol,:);
            
            xLabs{cueSide} = cueLabs{cueRow};
            
        case 2 %focal left
            cueRow = 2;
            sideCol = 1;
            
            allEffects(cueSide,1,:) = targEffects(cueRow, sideCol, :); %left word effect
            allEffects(cueSide,2,:) = distEffects(cueRow, sideCol, :); %right word effect
            
            xLabs{cueSide} = [cueLabs{cueRow} ' ' sideLabs{sideCol}];
            
        case 3 %focal right
            cueRow = 2;
            sideCol = 2;
            
            allEffects(cueSide,1,:) = distEffects(cueRow, sideCol, :); %left word effect
            allEffects(cueSide,2,:) = targEffects(cueRow, sideCol, :); %right word effect
            
            xLabs{cueSide} = [cueLabs{cueRow} ' ' sideLabs{sideCol}];
    end
end

%% run stats on effects!

fprintf(1,'\n---------------------------------\n');
fprintf(1,'ANALYSIS OF EFFECTS OF FREQUENCY (HIGH-LOW) OF EACH WORD IN EACH CUE CONDITION\n');
fprintf(1,'In focal cue conditions, the effect of the uncued word is from the "uncuedDprime" variable\n');

%make a table
nRows = numel(allEffects);

effects = NaN(nRows,1);
cueConds = cell(nRows,1);
wordSides = cell(nRows,1);
subjects = NaN(nRows,1);

effectCIs = NaN(3,2,2);
effectBootPs = NaN(3,2);
rowI = 0;
for cueI = 1:3
    for sideI = 1:2
        for si = 1:rAvg.nsubj
            rowI = rowI+1;
            effects(rowI) = allEffects(cueI, sideI, si);
            cueConds{rowI} = xLabs{cueI};
            wordSides{rowI} = sideLabs{sideI};
            subjects(rowI) = si;
        end
        %test this effect
        theseEs = squeeze(allEffects(cueI, sideI, :));
        [eCI, ~, eDist] = boyntonBootstrap(@nanmean,theseEs,1000,95,true);
        effectCIs(cueI, sideI,:) = eCI;
        effectBootPs(cueI, sideI) = getBootPs(eDist',0);
        
        [~,tp,~,tst] = ttest(theseEs);
        
        fprintf(1,'\nEffect of %s word in %s condition:\n', sideLabs{sideI}, xLabs{cueI});
        fprintf(1,'\tmean=%.4f, SEM = %.3f, 95%%CI = [%.3f %.3f], bootP=%.5f,', nanmean(theseEs), standardError(theseEs'), eCI(1), eCI(2), effectBootPs(cueI, sideI))
        fprintf(1,'\n\tt(%i)=%.3f,p=%.5f\n', tst.df, tst.tstat, tp);
        
    end
end

T2 = table;
T2.effect = effects;
T2.cue = cueConds;
T2.wordSide = wordSides;
T2.subject = categorical(subjects);

fprintf(1,'\n\nLME analysis of frequency effects, as a function of word side and cue condition:\n');

eqtn = 'effect ~ cue*wordSide';
lme = fitlme(T2, eqtn);
display(lme);
display(lme.anova);

%% Plot effects of each word (left, right) in each cue condition and
msToPlot = nanmean(allEffects,ndims(allEffects));
esToPlot = standardError(allEffects, ndims(allEffects));

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

opt.ylims = [-1 1];
opt.yLab = '\Delta d''';

opt.xTickLabs = {}; %xLabs';

opt.legendLabs = {'Left high-low', 'Right high-low'};
opt.legendLoc = 'SouthEast';
opt.lev1ForLegend = 2;


figure(figH);
subplot(1,2,1);

barCenters = barPlot_AW(msToPlot, esToPlot, opt);
%add stars for significance
for cueI = 1:3
    for sideI=1:2
        pval = effectBootPs(cueI, sideI);
        if pval<0.05
            starY = opt.ylims(1)+0.91*diff(opt.ylims);
            if pval<0.001
                startTxt = '***';
            elseif pval<0.01
                startTxt = '**';
            elseif pval<0.05
                startTxt = '*';
            end
            text(barCenters(cueI,sideI),starY,startTxt,'HorizontalAlignment','center');
        end
    end
end

title('Task accuracy');
set(gca,'TitleFontSizeMultiplier',1,'TitleFontWeight','normal');



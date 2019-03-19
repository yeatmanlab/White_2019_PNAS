%% Figure 5 for White, Palmer, Boynton & Yeatman, PNAS 2019
% Plots effects of lexical frequency on d' and on mean BOLD in left VWFA-2.
% 
% This function calls two other functions to make panels 5A and 5B. 
% by Alex L. White at the University of Washington, 2019


%% Figure 5: Effects of lexical frequency 
function Fig5_LexicalFrequency()

figH = figure; 

Fig5A_DprimeSidexFrequency(figH);
Fig5B_LVWFA2_FreqEffects(figH)

p = getPaths(); 

figSize = [9 4];
fontSize = 9;

set(gcf,'units','centimeters','pos',[8 8 figSize],'color','w');
figTitle = fullfile(p.figures,sprintf('Fig5_FrequencyEffects.eps'));
exportfig(gcf,figTitle,'Format','eps','bounds','loose','color','rgb','LockAxes',0,'FontMode','fixed','FontSize',fontSize);

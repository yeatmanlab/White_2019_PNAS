%% function [pProcessBoth, pTask1First, slope, intercept, distFromAllOrNone] = AnalyticDualTaskSerialModel(singleAccs, dualAccs, doPlot)
% by Alex White, 2017
%
% Given accuracy (p(correct) or area under ROC curve) in dual and single-task 
% conditions for two tasks (e.g., left and right) this model finds the best-fitting
% generalized serial model. This is an analytic solution, by which I mean
% that it does not do any optimization search or simulate individual
% trials. It basically just re-parameterizes the dual-task accuracy levels
% into two other parameters, given single-task accuracy levlels. 
% 
% The model assumes that on each trial one task (or side) is attended 'first', and the 
% other second. On some fraction of trials, both get processed, and on the 
% remaining trials, only the first one is processed and the subject must guess about the other. 
% The proportion of trials when both are processed is called pProcessBoth.
% The proportion of trials when task number 1 is attended first is called
% pTask1First.
% 
%
% This function returns those two free parameters as well as a slope and
% intercept for the best-fitting line in the AOC plot. That line is
% determined completely by the data and the processBoth parameter. 
% 
% Inputs: 
% - singleAccs: a 1x2 vector of accuracy in single-task conditions for the
%   two tasks. Units: p(correct) or area under ROC curve. 
% - dualACcs: a 1x2 vector of accuracy in dual-task conditions for the same
%   two tasks.
% - doPlot: whether to make an ROC plot showing data and model predictions.
% 
% Outputs: 
% - pProcessBoth: the proportion of trials when the subject successfully processed both 
%   stimuli [each with   p(correct) = single-task p(correct)]. On the remaining 
%   trials, the subject is forced completely guess about the task not attended 'first'. 
% - pTask1First: best-fitting parameter for proportion of dual-task trials when task number 1
%   is attended 'first'. 
% - slope: slope of the best-fitting line in the AOC plot
% - intercept: y-intercept of the best-fitting line in the AOC plot 
% - distFromAllOrNone: distance between dual-task data point and the
%   nearest point on the all-or-none switching model. Negative if data point
%   is below the line. 
%
function [pProcessBoth, pTask1First, slope, intercept, distFromAllOrNone] = AnalyticDualTaskSerialModel(singleAccs, dualAccs, doPlot)

%We assume that chance level is 0.5, and the axes of the AOC plot are
%limited to [0.5 1]. So to calculate things like slopes and intercepts and
%distances in the AOC space, we need to subtract 0.5 from everything and
%pretend like the axis limits really are [0 0.5]. 
singleAccs = singleAccs - 0.5; 
dualAccs = dualAccs - 0.5; 


%pull out data:
dual1  = dualAccs(1); %y-value
dual2 = dualAccs(2); %x-value
single1  = singleAccs(1);%y-value
single2 = singleAccs(2);%x-value


%Given single-task accuracy on the two tasks, single1 and single2,
%and dual-task accuries, dual1 and dual2,
%we want to solve for the dual-task model parameters pProcessBoth and
%pTask1First.

%First we can solve for the slope of all the dual-task lines
%slope = (.5-single1)/(single2-0.5); %this was true before we subtracted 0.5 from all data points 
slope = -single1/single2;

%Second, we can solve for the y-intercept of the best-fitting line, using our measured
%left and right dual-task performance
intercept = dual1 - slope*dual2;

%Then we can find the predicted value of left-dual task performance when attention is
%devoted entirely to task 2, and task 1 is supposed to be ignored on 100%
%of trials. Call this dual1Ignored (dual-task one when ignored). It is the y-value of the line when the x-value is
%single2.
dual1Ignored = slope*single2+intercept;

% from this we can solve for the proportion of trials when only 1 side is
% processed and the other must be guessed. 
% We assume that accuracy when guessing is 0.5.
guessingPCorr = 0; %0.5; %no longer 0.5 because we subtracted that out 
pProcessOnlyOne = (dual1Ignored - single1)/(guessingPCorr - single1);
%clip this propability to be between 0 and 1
pProcessOnlyOne(pProcessOnlyOne<0)=0; 
pProcessOnlyOne(pProcessOnlyOne>1)=1; 

pProcessBoth  = 1-pProcessOnlyOne;

%[We could also do this using the predicted right-side accuracy when right side is
%supposed to be ignored, dual2Ignored
%dual2Ignored = (single1-intercept)/slope;
%pProcessOnlyOne_B = (dual2Ignored - single2)/(guessingPCorr - single2); %]

% Then solve for pAttendTask1, or pAL (proportion of dual-task trials
% attend to left first)
pTask1First = (dual2 - single2)/(pProcessOnlyOne*guessingPCorr - single2*pProcessOnlyOne);

%clip in case dual-task performance exceeds any serial model 
pTask1First(pTask1First>1) = 1;
pTask1First(pTask1First<0) = 0; 


%% find the minimum distance between the dual-task data point and the all-or-none serial prediction line 
%difference in y-intercepts of best-fitting line and all-or-none line
% dY = intercept - single1;
% 
% %compute x-intercept of best-fitting line
% xIntercept = -1*intercept/slope;
% %difference between this and the all-or-none's x-intercept (which is simple
% %single-task accuracy level on x-axis)
% dX = xIntercept - single2;
% 
% %hypotenuese
% h = sqrt(dX^2 + dY^2); 
% %some trig:
% distFromAllOrNone1 = sqrt(dX^2 - (h/2)^2);

distFromAllOrNone = distanceFromPointToLine(dual2,dual1,slope,single1);


%if performance is WORSE than all-or-none model, make this negative
%difference in y-intercepts of best-fitting line and all-or-none line
dY = intercept - single1;
if dY<0
    distFromAllOrNone = distFromAllOrNone*-1;
end


%% plot
if doPlot
    dataMarkSz = 30;
    modelMarkSz = 10;
    
    figure; hold on;
    
    %plot box constrained by indepenent processing:
    plot([single2 single2],[0.5 single1],'k-');
    plot([0.5 single2],[single1 single1],'k-');
    
    %plot all-or-none serial model. Solve for that equation using x,y
    %coordinates of left single task accuracy:
    x0 = 0.5; y0=single1;
    betaAllOrNone = y0 - slope*0.5;
    xs = [.5 single2];
    ys = slope*xs+betaAllOrNone;
    plot(xs,ys,'k-');
    
    %plot the best-fitting serial model line (fit with pProcessOnlyOne, still free
    %parameter of pAttendTask1)
    dual2Ignored = (single1-intercept)/slope;

    xs = [dual2Ignored single2];
    ys = slope*xs+intercept;
    plot(xs,ys,'r-');
    
    
    %Plot prediction of fixed-capacity parallel processing
    singleAcc = [single2 single1];
    singleDs = sqrt(2)*norminv(singleAcc); %convert to d' (assuming A' is like PC?')
    pSamplesOnTask1 = 0:0.01:1;  %vector of proportion of fixed number of sensory "samples" devited to task 1
    
    % How to calculate d' for diffent numbers of "samples"
    % given sample sizes ss1 and ss2: ss1 = ss2/2
    % std1 = sqrt(2)*std2
    % generalized:
    %if ss1 = ss2*q, then
    %std1 = sqrt(1/q)*std2
    
    %therefore, d' relation is:
    % d1 = d2/sqrt(1/q)
    
    %dual-task dprimes
    fixedD1s = singleDs(1)./sqrt(1./pSamplesOnTask1);
    fixedD2s = singleDs(2)./sqrt(1./(1-pSamplesOnTask1));
    
    %dual-task AROCs,
    fixedA1s = normcdf(fixedD1s/sqrt(2));
    fixedA2s = normcdf(fixedD2s/sqrt(2));
    
    plot(fixedA1s,fixedA2s,'y-','LineWidth',2)
    
    
    %plot data:
    plot(0.5,single1,'b.-','MarkerSize',dataMarkSz);
    plot(single2,0.5,'b.-','MarkerSize',dataMarkSz);
    plot(dual2,dual1,'k.','MarkerSize',dataMarkSz)
    
    
    %plot the predicted dual-task performance given our 2 parameters:
    predDual = SerialDualTaskAccGivenSingleTask([single1 single2],pTask1First,pProcessBoth);
    plot(predDual(2), predDual(1),'r.','MarkerSize',modelMarkSz);

    axis square;
    xlabel('Right side p(c)');
    ylabel('Left side p(c)');
    xlim([0.5 1]); ylim([0.5 1]);
    ticks = 0.5:0.1:1;
    set(gca,'XTick',ticks,'YTick',ticks);
    
end
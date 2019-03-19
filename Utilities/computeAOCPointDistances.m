%function computeAOCPointDistances(as,es,plotOpt)
%
%Inputs:
%- as: 2x2 matrix of BOLD response differences, relative to the 'ignored' responses.
%    rows = single-task, dual-task
%    columns = task 1 (y-axis), task 2 (x-axis)
%- pTask1First: estimated from behavioral accuracy data, this is the
% probability that the subject processed task 1 first, assuming a generalized serial model is true
%
%
% Outputs:
% - distFromAllOrNoneLine: distance of the dual-task data point from the serial switching model line.
%    Positive values mean above the line, negative below.
% - distFromAllOrNoneAccPt: distance of the dual-task data point from the
% point along the swerial switching line predicted by  pTask1First. That
% predicted point is simply pTask1First proportion of the way between the
% two single-task points on the axes. Postive values are for points above
% the serial switching line, negative for points below.
% - distFromUnlimited: distance of the dual-task point from the
%   unlimited capacity independence point. Negative if it's below the line with slope = -1
%   that passes through the unlimited capacity point. %hat's true if the average of
%  (dual1-single1) and (dual2-single2) is negative.
%
% by Alex White, 2018

function [distFromAllOrNoneLine, distFromAllOrNoneAccPt, distFromUnlimited] = computeAOCPointDistances(as,pTask1First)


%% estimate distance from all-or-none

%pull out data:
dualAccs = as(2,:);
singleAccs = as(1,:);

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


distFromAllOrNoneLine = distanceFromPointToLine(dual2,dual1,slope,single1);


%if performance is WORSE than all-or-none model, make this negative
%difference in y-intercepts of best-fitting line and all-or-none line
dY = intercept - single1;
if dY<0
    distFromAllOrNoneLine = distFromAllOrNoneLine*-1;
end

%% estimate distance from unlimited capacity point

distFromUnlimited = sqrt((single1-dual1)^2 + (single2-dual2)^2 );

%make it negative if it's below the line with slope = -1 that passes
%through the unlimited capacity point.
%that's true if the average of (dual1-single1) and (dual2-single2) is
%negative.

if mean([(dual1-single1) (dual2-single2)])<0
    distFromUnlimited = -1*distFromUnlimited;
end

%% esetimate distance from the point predicted by the subject's behavioral accuracy -
%that is, the estimated probability of procesesing stimulus 1 first in the
%dual-task condition

if ~isempty(pTask1First)
    
    predPt = [0 single1].*pTask1First + [single2 0].*(1-pTask1First);
    
    distFromAllOrNoneAccPt = sqrt((dual1-predPt(2))^2 + (dual2-predPt(1))^2);
    
    %if performance is WORSE than all-or-none model, make this negative
    %difference in y-intercepts of best-fitting line and all-or-none line
    dY = intercept - single1;
    if dY<0
        distFromAllOrNoneAccPt = distFromAllOrNoneAccPt*-1;
    end
    
else
    distFromAllOrNoneAccPt = NaN;
end



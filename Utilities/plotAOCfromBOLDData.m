%function plotAOCfromBOLDData(as,es,plotOpt)
%
%Inputs:
%- as: 2x2 matrix of BOLD response differences, relative to the 'ignored' responses.
%    rows = single-task, dual-task
%    columns = task 1 (y-axis), task 2 (x-axis)
%- es: a 2x2 matrix of error bars for each of the conditions in as
%- plotOpt: structure with fields including::
%      - edgeColors: a 2x3 matrix with each row a color for marker edge.
%                    rows = single-task, dual-task
%      - fillColors: like edgeColors, but for fill of markers
%      - markSz: size of markers
%      - doLegend: whether to make legend
%      - doXLabel: whether to add x-axis label
%      - doYLabel: whether to add y-axis label
%      - sideLabels: a 1x2 cell array for the labels of the two axes
%
% Outputs: 
% - distFromAllOrNone: distance of the dual-task data point from the serial switching model line.  
%    Positive values mean above the line, negative below 
% - distFromUnlimited: distance of the dual-task point from the
%   unlimited capacity independence point. Negative if it's below the line with slope = -1 
%   that passes through the unlimited capacity point. That's true if the average of 
%  (dual1-single1) and (dual2-single2) is negative. Positive otherwise. 
%
%
% by Alex White, 2018

function [distFromAllOrNone, distFromUnlimited] = plotAOCfromBOLDData(as,es,plotOpt)

if ~isfield(plotOpt,'markSz')
    plotOpt.markSz = 12;
end
if ~isfield(plotOpt,'doLegend')
    plotOpt.doLegend = true;
end
if ~isfield(plotOpt,'doXLabel')
    plotOpt.doXLabel = true;
end
if ~isfield(plotOpt,'doYLabel')
    plotOpt.doYLabel = true;
end

if ~isfield(plotOpt,'errorBarWidth')
    plotOpt.errorBarWidth = 2;
end

if ~isfield(plotOpt,'sideLabels')
    plotOpt.sideLabels = {'Right hem. dPSC','Left hem. dPSC'};
end

singleTaskGroundLevel = 0; %for accuracy data, this would be 0.5

%% setup axes
if ~isfield(plotOpt,'axlims')
    axlims = [0.5 1];
else
    axlims = plotOpt.axlims;
end

if ~isfield(plotOpt,'axticks')
    axticks = linspace(axlims(1),axlims(2),5);
else
    axticks = plotOpt.axticks;
end
if ~isfield(plotOpt,'axtickLabels')
    axtickLabels = cell(1,length(axticks));
    for ai=1:length(axticks)
        v = axticks(ai);
        if mod(ai,2)==1  %if ai==1 || ai==length(axticks)
            if (round(v)-v)==0
                frmt = '%.1f';
            elseif (round(10*v)-10*v) == 0
                frmt = '%.1f';
            elseif (round(100*v)-100*v) == 0
                frmt = '%.2f';
            else
                frmt = '%.3f';
            end
            axtickLabels{ai} = sprintf(frmt,v);
        else
            axtickLabels{ai} = '';
        end
    end
else
    axtickLabels = plotOpt.axtickLabels;
end
markSz = plotOpt.markSz;
axLineWidth = plotOpt.axLineWidth;
datLineWidth = plotOpt.datLineWidth;
errorBarWidth = plotOpt.errorBarWidth;

%% set colors
singEdgeColor = plotOpt.edgeColors(1,:);
singFillColor = plotOpt.fillColors(1,:);
dualEdgeColor = plotOpt.edgeColors(2,:);
dualFillColor = plotOpt.fillColors(2,:);

%% plot 0 lines 
hold on;

if any(axlims<0)
    plot([0 0],axlims,'k-'); 
    plot(axlims,[0 0],'k-');
end
%% Plot the box for unlimited capacity model

plot([as(1,2) as(1,2)],[singleTaskGroundLevel as(1,1)],'k--','LineWidth',datLineWidth);
plot([singleTaskGroundLevel as(1,2)], [as(1,1) as(1,1)],'k--','LineWidth',datLineWidth);

%% Plot the serial model prediction of straight line:
plot([singleTaskGroundLevel as(1,2)], [as(1,1) singleTaskGroundLevel],'k-','LineWidth',datLineWidth);


%% Plot data
%First task single-task performance (side==1) on y axis
%add error bar:
plot(singleTaskGroundLevel([1 1]),as(1,1)+[-1 1]*es(1,1),'-','Color',singEdgeColor,'LineWidth',errorBarWidth);
hSing=plot(singleTaskGroundLevel(1),as(1,1),'o','MarkerSize',markSz,'MarkerFaceColor',singFillColor,'MarkerEdgeColor',singEdgeColor);

%Second task single-task performance (side==2) on x axis
%add error bar:
plot(as(1,2)+[-1 1]*es(1,2),singleTaskGroundLevel([1 1]),'-','Color',singEdgeColor,'LineWidth',errorBarWidth);
plot(as(1,2),singleTaskGroundLevel(1),'o','MarkerSize',markSz,'MarkerFaceColor',singFillColor,'MarkerEdgeColor',singEdgeColor);

%dual-task performance
%with error bars:
dualx = as(2,2); dualy=as(2,1);
dualxE = es(2,2); dualyE = es(2,1);
plot(ones(1,2)*dualx,dualy+[-1 1]*dualyE,'-','Color',dualEdgeColor,'LineWidth',errorBarWidth);
plot(dualx+[-1 1]*dualxE,ones(1,2)*dualy,'-','Color',dualEdgeColor,'LineWidth',errorBarWidth);

hDual = plot(dualx,dualy,'o','MarkerSize',markSz,'MarkerFaceColor',dualFillColor,'MarkerEdgeColor',dualEdgeColor,'LineWidth',datLineWidth+2);



xlim(axlims);
ylim(axlims);
set(gca,'XTick',axticks,'XTickLabel',axtickLabels);
set(gca,'YTick',axticks,'YTickLabel',axtickLabels);
set(gca,'LineWidth',axLineWidth);
set(gca,'FontName','Helvetica')

if plotOpt.doXLabel
    xlabel(plotOpt.sideLabels{2});
end
if plotOpt.doYLabel
    ylabel(plotOpt.sideLabels{1});
end
axis square;

if plotOpt.doLegend
    legend([hSing hDual],{'Focal cued - focal uncued','Distributed - focal uncued'},'Location','NorthWest'); legend boxoff;
end

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

distFromAllOrNone = distanceFromPointToLine(dual2,dual1,slope,single1);

%if performance is WORSE than all-or-none model, make this negative
%difference in y-intercepts of best-fitting line and all-or-none line
dY = intercept - single1;
if dY<0
    distFromAllOrNone = distFromAllOrNone*-1;
end

%% estimate distance from unlimited capacity point 

distFromUnlimited = sqrt((single1-dual1)^2 + (single2-dual2)^2);

%make it negative if it's below the line with slope = -1 that passes
%through the unlimited capacity point. 
%that's true if the average of (dual1-single1) and (dual2-single2) is
%negative. 

if mean([(dual1-single1) (dual2-single2)])<0
    distFromUnlimited = -1*distFromUnlimited;
end



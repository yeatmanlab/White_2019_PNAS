%% function dualPCorr = SerialDualTaskAccGivenSingleTask(singlePCorr, pTask1First, pProcessBoth) 
% 
% Implements a serial model of dual-task performance.  This is an analyitic solution to
% predict dual-task accuracy given single-task performance (rather than
% simulating tons of trials).
% Specifically, it computes p(correct) for two tasks, given single-task accuracies,
% proportion of dual-task trials that task 1 is done 'first', and the
% proportion of trials when both tasks can be completed. 
%  
% This is a serial model, meaning that only 1 task is done at a time. 
% We assume that on each trial, 1 task is attended first. That task will be completed, yielding the same
% p(corr) as in the single-task condition. The ignored task will either not
% be processed at all, in which case the subject must guess about it, or fully processed. 
% The probability of getting to do the second task is pProcessBoth.

% 
% Inputs: 
% - singlePCorr: a 1x2 vector of p(corr) for the single-task conditions. 
% - pTask1First: the proportion of dual-task trials that task 1 is
%   'attended' 'first'
% - pProcessBoth: the proportion of trials when the subject successfully processed both 
%   stimuli. On the remainder, the subject is forced completely guess about
%   the task not attended 'first'. 
% 
% by Alex White, 2017


function dualPCorr = SerialDualTaskAccGivenSingleTask(singlePCorr, pTask1First, pProcessBoth) 

pStartWithTaskI = [pTask1First (1-pTask1First)];

guessingPCorr = 0.5;

dualPCorr = NaN(1,2); 
for i=1:2
    pProcess = pProcessBoth + (1-pProcessBoth)*pStartWithTaskI(i);
    dualPCorr(i) = pProcess*singlePCorr(i) + (1-pProcess)*guessingPCorr;
end

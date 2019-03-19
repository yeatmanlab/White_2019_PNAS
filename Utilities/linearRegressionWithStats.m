%% function [betas, rSqr, r2adj, SSres, P] = linearRegressionWithStats(X, Y) 
% by Alex White, 2018
% Solve the linear model: Y = X*betas, given observered X and Y. 
% Returns the least-squares solution for betas, and also r-squared,
% adjusted r-squared, the sum of squared errors of the model fit, and the
% number of free parametres. 
% Inputs: 
% - X: a n x c matrix of predictor values 
% - Y: a n x j matrix of observed values 
% 
% Outputs: 
% - B: a c x j matrix of beta weights 
% - rSqr: the r-squared value for the fit: 1-SSres/SStot
% - rSqrAdj: r-squared adjusted bye the degrees of freedom: 
%   1-(1-rSqr)*(n-1)/(n-P-1) 
%   where n is the number of data points and P is the number of free
%   parameters
% - SSres: sum of squared errors between predicted values of Y and actual Y
% - P: number of free parameters (numel(B)). 

function [B, rSqr, rSqrAdj, SSres, P] = linearRegressionWithStats(X, Y) 

if size(X,1) ~= size(Y,1) 
    error('X and Y must have the same number of rows');
end

B = pinv(X)*Y;

%Compute r-squared
Yhat = X*B;
resids = Y - Yhat;
meanB = nanmean(Y(:));
SSres = sum(resids(:).^2);
SStot = sum((Y(:) - meanB).^2);
rSqr = 1-SSres/SStot;

%compute adjusted R2:
N = size(X,1); %number of data points
P = numel(B); %number of free parameters

if N>(P+1) %only compute if there are more data points than free parameters
    rSqrAdj = 1 -(1-rSqr)*(N-1)/(N-P-1);
else
    rSqrAdj = NaN;
end

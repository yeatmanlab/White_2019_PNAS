% function d = distanceFromPointToLine(x,y,m,k)
% from https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Another_formula
% 
% Inputs: 
% x,y are coordinates of point 
% m = slope of line 
% k = y-intercept of line 
% 
% Outputs: 
% - d = distance between point and line 

function d = distanceFromPointToLine(x,y,m,k)
if isinf(m) || isnan(m)
    fprintf(1,'\n(distanceFromPointToLine) Warning: slope m is inf or undefined. Setting distance to NaN\n');
    d= NaN;
else
    d = abs(k + m*x-y) / sqrt(1+m^2); 
end



function [advecValue, beta] = interpolateForValue(localInterpPoints, M, newLocalPointX, values)
%% INTERPOLATEFORVALUE() uses least squares to interpolate through scalar values
%
% interpolateForValue uses quadratic least squares to fit a parabola
% through a set of scalar values associated with the surface. The function
% also ensures that that new value is within the interpolation bounds
%
% INPUTS
%   localInterpPoints       ===     struct array: contains current point neighbours, see neighbourSearch function header for properties
%   M                       ===     matrix(3x3): constructed of x-values used for interpolation
%   newLocalPointX          ===     scalar: x coordinate of newly resampled point in local coordinates
%
%
% Author: Solene Hegarty-Cremer
%%

% %Use values as y points
b = [sum(values);...
    sum([localInterpPoints.xCoord].*values);...
    sum(([localInterpPoints.xCoord].^2).*values)];

beta = M\b;

%Solves advection
advecValue = beta(1) + beta(2)*newLocalPointX + beta(3)*newLocalPointX^2;

if advecValue > max(values)
    advecValue = max(values);
elseif advecValue < min(values)
    advecValue = min(values);
end

end

function [a, M, b] = quadInterpolate(localInterpPoints)
%% QUADINTERPOLATE() use least squares to interpolate through local points
%
% quadInterpolate uses quadratic least squares to fit a parabola through
% the local points. This results in a local representation of the position
% of the surface.
%
% INPUTS
%   localInterpPoints    ===     struct array: contains current point neighbours, see neighbourSearch function header for properties
%
%
% Author: Solene Hegarty-Cremer
%%
xSum = sum([localInterpPoints.xCoord]);
xSumSquare =  sum([localInterpPoints.xCoord].^2);
xSumCube =  sum([localInterpPoints.xCoord].^3);

M = [length(localInterpPoints), xSum, xSumSquare;...
    xSum, xSumSquare, xSumCube;...
    xSumSquare, xSumCube, sum([localInterpPoints.xCoord].^4)];

b = [sum([localInterpPoints.yCoord]);...
    sum([localInterpPoints.xCoord].*[localInterpPoints.yCoord]);...
    sum(([localInterpPoints.xCoord].^2).*[localInterpPoints.yCoord])];

a = M\b;
end
function [interpPoints, toRemove] = findMNeighbours(currentPoint, activePoints, m, x, y, domainLen, thetaMin, dx)
%% FINDMNEIGHBOURS() finds the closest m marker particles to the current point
%
% findMNeighbours conducts a local search along the underlying Eulerian
% grid to find the m+2 closest neighbouring marker particles to the current
% point of interest and stores their values of interest in a new struct.
% The points are ordered according to their distance from the current
% marker particle and points are disregarded if their normals vary too
% greatly from the normal of the current point.
%
% INPUTS
%   currentPoint    ===     struct: single struct object of point being evaluated, see CBPM function header for properties
%   activePoints    ===     struct array: contains all domain grid cells, see CBPM function header for properties
%   x               ===     vector(1xdomainLen): x locations of grid cell edges
%   y               ===     vector(1xdomainLen): y locations of grid cell edges
%   domainLen       ===     scalar: length of one direction of discretised domain
%   m               ===     vector(1x2): number of neighbours to use for interpolation for surface position and scalar quantities (density) respectively
%   thetaMin        ===     scalar: parameter used to detect collisions
%   dx              ===     scalar: space discretisation (the same in x and y directions)
%
%
% Author: Solene Hegarty-Cremer
%%
% Initialise interpPoints
j = 1;
interpPoints(1:m+2) = struct('gridPointCoords', [],...
    'footPointCoords', [], ...
    'normal', [],...
    'val', [], ...
    'vs', [],...
    'cellType', []);
%Current point is first closest point
interpPoints(j) = struct('gridPointCoords', currentPoint.gridPointCoords,...
    'footPointCoords', currentPoint.footPointCoords, ...
    'normal', currentPoint.normal,...
    'val', currentPoint.val, ...
    'vs', currentPoint.vs,...
    'cellType', currentPoint.cellType);

%If footpointcoords(1) == 100 this is a dummy point from activation (don't
%include in neighbours)
if interpPoints(1).footPointCoords(1) ~= 100
    j = 1; %number of found neighbours
else
    j = 0;
end

%Set window already searched
xLeft = currentPoint.gridPointIndices(1);
xRight = currentPoint.gridPointIndices(1);
yDown = currentPoint.gridPointIndices(2);
yUp = currentPoint.gridPointIndices(2);

%Find at least m neighbours
[interpPoints, xLeft, xRight, yDown, yUp, toRemove]  = neighbourSearch(interpPoints, activePoints, ...
    domainLen, x, y, dx, xLeft, xRight, yDown, yUp, m, j);

%Get distances
distances = zeros(length([interpPoints.val]),1);
for k = 1:length(distances)
    %Distance away from current grid point
    distances(k) = norm(currentPoint.gridPointCoords - interpPoints(k).footPointCoords);
end

%Sort yi in ascending order according to distance
[~, perm] = sort(distances);
interpPoints = interpPoints(perm);

%Check normals
pointsToRemove = [];
for k = 2:length(interpPoints)
   if dot(interpPoints(1).normal, interpPoints(k).normal) < cos(thetaMin)
       pointsToRemove = [pointsToRemove, k];
   end
end
interpPoints(pointsToRemove) = [];

%If we now have less than m neighbours search again
while length(interpPoints) < m && toRemove ~= 1
    pointsToRemove = [];
    [interpPoints, xLeft, xRight, yDown, yUp, toRemove]  = neighbourSearch(interpPoints, activePoints, ...
        domainLen, x, y, xLeft, xRight, yDown, yUp, m, length(interpPoints));
    
    for k = 2:length(interpPoints)
        if dot(interpPoints(1).normal, interpPoints(k).normal) < cos(thetaMin)
            pointsToRemove = [pointsToRemove, k];
        end
    end
    interpPoints(pointsToRemove) = [];
    
end

if toRemove ~= 1
    interpPoints = interpPoints(1:m);
end
end
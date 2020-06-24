function [activePoints, newPoints] = activateNewPoints(activePoints, m, x, y, dx, domainLen, thetaMin, Vx, Vy, t)
%% ACTIVATENEWPOINTS actives grid cells which newly contain a portion of the interface
%
% activateNewPoints considers inactive grid cells neighbouring active ones
% for activation. The direction the interface is moving is considered to
% narrow the list of inactive cells to consider. If these cells now contain
% a portion of interface, they are activated with marker particles.
%
% INPUTS
%   activePoints    ===     struct array: contains all domain grid cells, see CBPM function header for properties
%   m               ===     vector(1x2): number of neighbours to use for interpolation for surface position and scalar quantities (density) respectively
%   x               ===     vector(1xdomainLen): x locations of grid cell edges
%   y               ===     vector(1xdomainLen): y locations of grid cell edges
%   domainLen       ===     scalar: length of one direction of discretised domain
%   thetaMin        ===     scalar: parameter used to detect collisions
%   Vx              ===     anon. function(point, t): describes change in x position of the interface, takes in point struct and time
%   Vy              ===     anon. function(point, t): describes change in y position of the interface, takes in point struct and time
%   t               ===     scalar: current simulation time
%
%
% Author: Solene Hegarty-Cremer
%%
%Count new points
newPoints = [];

%Gather points to consider
pointsToCheck = getInactiveNeighbours(activePoints, domainLen, x, y, Vx, Vy, t);

%For each of these, calculate grid point
for i = 1:size(pointsToCheck, 1)
    currentPoint = struct('gridPointCoords', [x(pointsToCheck(i,1)), y(pointsToCheck(i,2))],...
        'gridPointIndices', pointsToCheck(i,:), ...
        'footPointCoords', [100 100], ...
        'normal', [1, 0], ...
        'val', 0, ...
        'curv', 0, ... 
        'vs', 0,...
        'cellType', 0);
   
    %Interpolate surface at proposed grid cell
    [newFootPoint, newNormal, newValue, curvature, ...
        dist, neighbours, curveLength, alpha, beta, localPoint, zeta, newV,...
        cellType, toRemove] = interpolateAndResampleOnePoint(currentPoint, activePoints, m,...
    x, y, domainLen, thetaMin, dx);
    
    %Activate if within grid cell and if point has not been flagged for
    %removal
    if  dist <= dx/2 && ~toRemove
        activePoints(pointsToCheck(i,1) + domainLen*(pointsToCheck(i,2)-1)) = struct('gridPointCoords', currentPoint.gridPointCoords,...
            'gridPointIndices', currentPoint.gridPointIndices, ...
            'footPointCoords', newFootPoint', 'dist', dist, 'normal', newNormal', ...
            'curv', curvature, 'val', newValue, 'newVal', newValue', ... 
            'newFootPointCoords', newFootPoint', 'neighbours', neighbours, 'curveLength', curveLength,...
            'alpha', alpha, 'beta', beta, 'localPoint', localPoint, 'zeta', zeta, 'vs', newV, ...
            'arcPoint', 0, 'cellType', cellType);
        
        %Keep track of newly activated points
        newPoints = [newPoints, pointsToCheck(i,1)+domainLen*(pointsToCheck(i,2)-1)];
    end
end
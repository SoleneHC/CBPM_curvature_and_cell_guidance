function [newFootPoint, newNormal, newValue, curvature, ...
    dist, neighbours, curveLength, alpha, beta, newLocalPoint, zeta, newV, ...
    newCellType, toRemove] = interpolateAndResampleOnePoint(currentPoint, ...
    activePoints, m,...
    x, y, domainLen, thetaMin, dx)
%% INTERPOLATEANDRESAMPLEONEPOINT Interpolates and resamples a single marker particle for interface position and scalar values
%
% interpolateAndResampleOnePoint uses local quadratic least square to
% interpolate for surface position and any scalar values associated with
% the marker particles (density, cellType, and vs). First, the closest m
% neighbours are found, then these are rotated into a local coordinate
% basis and the interpolation takes place. Next the resampling occurs. In
% any of these stages the point may be flagged for deactivation, finally
% the point may be deactivated if the curvature at that point exceeds
% maxCurvature or if the resampled location is outside of interpolation
% bounds.
%
% INPUTS
%   currentPoint    ===     struct: single struct object of point being evaluated, see CBPM function header for properties
%   activePoints    ===     struct array: contains all domain grid cells, see CBPM function header for properties
%   dx              ===     scalar: space discretisation (the same in x and y directions)
%   x               ===     vector(1xdomainLen): x locations of grid cell edges
%   y               ===     vector(1xdomainLen): y locations of grid cell edges
%   domainLen       ===     scalar: length of one direction of discretised domain
%   m               ===     vector(1x2): number of neighbours to use for interpolation for surface position and scalar quantities (density) respectively
%   thetaMin        ===     scalar: parameter used to detect collisions
%
%
% Author: Solene Hegarty-Cremer
%%
% curvature bound
maxCurvature = .25/dx;

% Find m neighbours
[interpPoints, toRemove] = findMNeighbours(currentPoint, activePoints, max(m), x, y, domainLen, thetaMin, dx);

if toRemove ~= 1
    % Change coordinates
    [localInterpPoints, localGridPoint, changeOfBasis] = changeCoordinates(interpPoints, currentPoint, dx);
    
    % Interpolate for interface position
    [alpha, M1, ~] = quadInterpolate(localInterpPoints(1:m(1)));
    
    % Resample for interface position
    [newLocalPoint, curvature, newLocalNormal, dist, neighbours, curveLength, toRemove] = localResample(alpha, ...
        localGridPoint);

    
    [~,M2,~] = quadInterpolate(localInterpPoints(1:m(2)));
    %Interpolate for value
    [newValue, beta] = interpolateForValue(localInterpPoints(1:m(2)), M2, newLocalPoint(1), [localInterpPoints(1:m(2)).val]);
    %Interpolate for vs
    [newV, zeta] = interpolateForValue(localInterpPoints(1:m(2)), M2, newLocalPoint(1), [localInterpPoints(1:m(2)).vs]);
    %Interpolate for cellType
    newCellType = interpolateForValue(localInterpPoints(1:m(1)), M1, newLocalPoint(1), [localInterpPoints(1:m(1)).cellType]);
    
    % Change back to Eulerian coords
    newFootPoint = changeOfBasis\(newLocalPoint) + interpPoints(1).footPointCoords';
    newNormal = changeOfBasis\newLocalNormal;
    dist = norm(newFootPoint' - currentPoint.gridPointCoords, 'inf');
    
    % Check deactivation conditions
    %Check for collision
    if collision(newFootPoint, newNormal, interpPoints, dx, thetaMin)
        toRemove = 1;
        fprintf('Point removed because of collision \n')
    end
    %Deactivate if curvature is too high or interpolation is out of
    %bounds
    if abs(curvature) > maxCurvature
        toRemove = 1;
        fprintf('Point disactivated because of curvature \n')
    elseif (newLocalPoint(1) < min(real([localInterpPoints.xCoord]))) || (newLocalPoint(1) > max(real([localInterpPoints.xCoord])))
        toRemove = 1;
        fprintf('Point disactivated because outside interpolation bounds \n')
    elseif toRemove == 1
        toRemove = 1;
    else
        toRemove = 0;
    end
else
    fprintf('Point removed because no neighbours')
    newFootPoint = 0;
    newNormal = 0;
    newValue = 0;
    curvature = 0;
    dist = 0;
    neighbours = 0;
    curveLength = 0;
    alpha = 0;
    beta = 0;
    zeta = 0;
    newV = 0;
    newCellType = 0;
    newLocalPoint = 0;
end
end
function [activePoints, removalPoints] = interpolateAndResample(activePoints, m, x, y, domainLen, thetaMin, dx)
%% INTERPOLATEAND RESAMPLE() interpolates and resamples the marker particles in activePoints
%
%
% interpolateAndResample applies the interpolation stage of the CBPM by
% interpolating and resampling each point in activePoints. It then flags
% the point to remove if needed otherwise updates all of its properties
%
%   activePoints    ===     struct array: contains all domain grid cells, see CBPM function header for properties
%   m               ===     vector(1x2): number of neighbours to use for interpolation for surface position and scalar quantities (density) respectively
%   x               ===     vector(1xdomainLen): x locations of grid cell edges
%   y               ===     vector(1xdomainLen): y locations of grid cell edges
%   domainLen       ===     scalar: length of one direction of discretised domain
%   thetaMin        ===     scalar: parameter used to detect collisions
%   dx              ===     scalar: space discretisation (the same in x and y directions)
%
% Author: Solene Hegarty-Cremer
%%
removalCount = 0;
removalPoints = [];

for i = find(~cellfun(@isempty, {activePoints.val})) %loop over active points   
        currentPoint = activePoints(i);
        
        [newFootPoint, newNormal, newValue, curvature, dist, ...
            neighbours, curveLength, alpha, beta, localPoint, zeta, newV,...
            newCellType, toRemove] = interpolateAndResampleOnePoint(currentPoint, ...
            activePoints, m, x, y, domainLen, thetaMin, dx);
        
        if toRemove
            removalCount = removalCount + 1;
            removalPoints(removalCount) = i;
        else
            activePoints(i).curv = curvature;
            activePoints(i).normal = newNormal';
            activePoints(i).dist = dist;
            activePoints(i).newFootPointCoords = real(newFootPoint');
            activePoints(i).newVal = newValue;
            activePoints(i).neighbours = neighbours;
            activePoints(i).curveLength = curveLength;
            activePoints(i).alpha = alpha;
            activePoints(i).beta = beta;
            activePoints(i).vs = newV;
%             activePoints(i).zeta = zeta;
            activePoints(i).cellType = newCellType;
            activePoints(i).localPoint = localPoint;
        end
end

% Update properties
for i = find(~cellfun(@isempty, {activePoints.val})) %over active points
    activePoints(i).footPointCoords = real(activePoints(i).newFootPointCoords);
    activePoints(i).val = activePoints(i).newVal;
end
end
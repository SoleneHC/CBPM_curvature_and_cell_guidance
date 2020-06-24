function activePoints = deactivatePoints(activePoints, removalPoints, dx)
%% DEACTIVATEPOINTS Deactivates points which have been previously flagged for removal
%
% deactivatePoints deactivates grid cells (removes their marker particles)
% for cells which have been flagged in this timestep for removal, also
% deactivates points which are too far away from their cell centers. This
% could be because the marker particle is no longer within the cell,
% because a collision was detected, or because the curvature at the given
% point exceeded the threshold.
%
% INPUTS
%   activePoints    ===     struct array: contains all domain grid cells, see CBPM function header for properties
%   removalPoints   ===     vector: indices of points flagged for removal
%
%
% Author: Solene Hegarty-Cremer
%%

activePoints(removalPoints) = struct('gridPointCoords', [], 'gridPointIndices', [], ...
                'footPointCoords', [], 'dist', [], 'normal', [], ...
                'curv', [], 'val', [], 'newVal', [], 'newFootPointCoords', [], ...
                'neighbours', [], 'curveLength', [],'alpha', [], 'beta', [],...
                'localPoint', [], 'vs', [], 'cellType', [],...
                'zeta', [], 'arcPoint', []);


for i = find(~cellfun(@isempty, {activePoints.val})) %over active points
    if activePoints(i).dist > dx/2 %norm(activePoints(i).gridPointCoords - activePoints(i).newFootPointCoords, 'Inf') > dx/1.8
        activePoints(i) = struct('gridPointCoords', [], 'gridPointIndices', [], ...
            'footPointCoords', [], 'dist', [], 'normal', [], ...
            'curv', [], 'val', [], 'newVal', [], 'newFootPointCoords', [], ...
            'neighbours', [], 'curveLength', [], 'alpha', [], 'beta', [],...
                'localPoint', [], 'vs', [], 'cellType', [],...
                'zeta', [], 'arcPoint', []);
    end
end
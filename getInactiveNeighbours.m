function pointsToCheck = getInactiveNeighbours(activePoints, domainLen, x, y, Vx, Vy, t)
%% GetInactiveNeighbours returns the set of indices of inactive cells neighbouring active ones
%
% getInactiveNeighbours gathers the indices of the set of inactive cells
% neighbouring active cells. The direction of the interface motion is used
% to construct this set (neighbours in the direction of flow are
% considered).
%
% INPUTS
%   activePoints    ===     struct array: contains all domain grid cells, see CBPM function header for properties
%   x               ===     vector(1xdomainLen): x locations of grid cell edges
%   y               ===     vector(1xdomainLen): y locations of grid cell edges
%   domainLen       ===     scalar: length of one direction of discretised domain
%   Vx              ===     anon. function(point, t): describes change in x position of the interface, takes in point struct and time
%   Vy              ===     anon. function(point, t): describes change in y position of the interface, takes in point struct and time
%   t               ===     scalar: current simulation time
%
%
% Author: Solene Hegarty-Cremer
%%
pointsToCheck = [];

for i = find(~cellfun(@isempty, {activePoints.val})) %over active points
    %If no neighbours add neighbour cross
    if size(activePoints(i).neighbours, 1) < 2
        activePoints(i).neighbours = [1 0; 0 1; -1 0; 0 -1];
        %Check if horizontal or vertical line - search more
    elseif all(activePoints(i).neighbours(1,:) == -activePoints(i).neighbours(2,:))
        %Add more to search
        activePoints(i).neighbours(3,:) = fliplr(activePoints(i).neighbours(1,:));
        activePoints(i).neighbours(4,:) = fliplr(activePoints(i).neighbours(2,:));
    end
    
    %Add diagonal neighbour in direction of flow
    activePoints(i).neighbours(end+1,:) = [-1*(Vx(activePoints(i), t) < 0) + 1*(Vx(activePoints(i),t) > 0), ...
        -1*(Vy(activePoints(i),t) < 0) + 1*(Vx(activePoints(i), t) > 0)];
    
    %Add diagonal neighbour in opposite direction of flow
    activePoints(i).neighbours(end+1,:) = [1*(Vx(activePoints(i), t) < 0) - 1*(Vx(activePoints(i),t) > 0), ...
        +1*(Vy(activePoints(i),t) < 0) - 1*(Vx(activePoints(i), t) > 0)];
    
    %Loop over cell neighbours
    for j = 1:size(activePoints(i).neighbours, 1)
        potentialNeighbour = activePoints(i).gridPointIndices + activePoints(i).neighbours(j,:);
        
        %Check in grid
        if potentialNeighbour(1) > 0 && potentialNeighbour(1) <= length(x) ...
                && potentialNeighbour(2) > 0 && potentialNeighbour(2) <= length(y)
            %Check inactive
            if isempty(activePoints(potentialNeighbour(1) + domainLen*(potentialNeighbour(2)-1)).val)
                pointsToCheck = [pointsToCheck; potentialNeighbour];
            end
            
        end
        
    end
    
    pointsToCheck = unique(pointsToCheck, 'rows');
end
function [interpPoints, xLeft, xRight, yDown, yUp, toRemove] = neighbourSearch(interpPoints, activePoints, ...
    domainLen, x, y, dx, xLeft, xRight, yDown, yUp, m, j)
%% NEIGHBOURSEARCH finds the m+2 closest neighbours which are active
%
% neighbourSearch searches outwards from the initial domain outlined by
% xLeft, xRight, yDown, yUp to find the m+2 closest active neighbours. As
% well as being active (having a marker particle inside the grid cell), the
% neighbours must also be far enough away from the current point to avoid
% any singularities in interpolation.
%
% INPUTS
%   interpPoints    ===     struct array: contains neighbours currently identified, properties:
%                               gridPointCoords
%                               footPointCoords
%                               normal
%                               val
%                               vs
%                               cellType
%   activePoints    ===     struct array: contains all domain grid cells, see CBPM function header for properties
%   x               ===     vector(1xdomainLen): x locations of grid cell edges
%   y               ===     vector(1xdomainLen): y locations of grid cell edges
%   dx              ===     scalar: space discretisation (the same in x and y directions)
%   domainLen       ===     scalar: length of one direction of discretised domain
%   m               ===     vector(1x2): number of neighbours to use for interpolation for surface position and scalar quantities (density) respectively
%   xLeft           ===     scalar: index of lower x boundary of domain already searched
%   xRight          ===     scalar: index of upper x boundary of domain already searched
%   yDown           ===     scalar: index of lower y boundary of domain already searched
%   yUp             ===     scalar: index of upper y boundary of domain already searched
%   j               ===     scalar: number of neighbours already found

%
% Author: Solene Hegarty-Cremer
%%
iters = 0;
maxIters = 3;
toRemove = 0;

while j < m+2 && iters <= maxIters
    % Establish points to search
    xPoints = max(xLeft-1, 1):min(xRight+1, length(x));
    yPoints = max(yDown-1, 1):min(yUp+1, length(y));
    
    topLine = [xPoints', ones(length(xPoints),1)*yPoints(1)];
    middle = [];
    for k = 2:length(yPoints)-1
        middle = [middle; xPoints(1), yPoints(k); xPoints(end), yPoints(k)];
    end
    bottomLine = [xPoints', ones(length(xPoints),1)*yPoints(end)];
    
    searchPoints = [topLine; middle; bottomLine];
    
    %Check if neighbours satisfy conditions
    for k = 1:size(searchPoints,1)
        if length(activePoints(searchPoints(k,1) + domainLen*(searchPoints(k,2)-1)).gridPointCoords) > 1  %if potential neighbour active
            if farEnoughAway(activePoints(searchPoints(k,1) + domainLen*(searchPoints(k,2)-1)), interpPoints, j, dx)
                j = j+1;
                interpPoints(j).gridPointCoords =  activePoints(searchPoints(k,1) + domainLen*(searchPoints(k,2)-1)).gridPointCoords;
                interpPoints(j).footPointCoords =  activePoints(searchPoints(k,1) + domainLen*(searchPoints(k,2)-1)).footPointCoords;
                interpPoints(j).normal =  activePoints(searchPoints(k,1) + domainLen*(searchPoints(k,2)-1)).normal;
                interpPoints(j).val = activePoints(searchPoints(k,1) + domainLen*(searchPoints(k,2)-1)).val;
                interpPoints(j).vs = activePoints(searchPoints(k,1) + domainLen*(searchPoints(k,2)-1)).vs;
                interpPoints(j).cellType = activePoints(searchPoints(k,1) + domainLen*(searchPoints(k,2)-1)).cellType;
            end
        end
    end
    
    xLeft = min(xPoints) - 1;
    xRight = max(xPoints) + 1;
    yDown = min(yPoints) - 1;
    yUp = max(yPoints) + 1;
    
    iters = iters + 1;
end

if iters > maxIters && j < m+2
    toRemove = 1;
end
end
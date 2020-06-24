function ok = farEnoughAway(potentialNeighbour, interpPoints, j, dx)
%% FARENOUGHAWAY() checks if potential neighbour is too close to interpolate
%
% INPUTS
%   interpPoints    ===     struct array: contains current point neighbours, see neighbourSearch function header for properties
%   j               ===     scalar: number of neighbours already found
%   dx              ===     scalar: space discretisation (the same in x and y directions)
%
%
% Author: Solene Hegarty-Cremer
%% 
ok = 1;

for i = 1:j
    if norm(potentialNeighbour.footPointCoords - interpPoints(i).footPointCoords, 'inf') < dx/2
        ok = 0;
    end
end


end
function [localInterpPoints, currentPointInfo, changeOfBasis] = changeCoordinates(interpPoints, currentPoint, dx)
%% CHANGECOORDINATES changes point in interp points to local coordinates
%
% changeCoordinates forms a local coordinate basis using the normal and
% tangent of the point closest to current point (interpPoints(1)). The rest
% of the interpolation points are rotated and translated to be in this
% coordinate basis. Also, the current point itself is also changed into the
% local coordinate basis.
%
% INPUTS
%   interpPoints    ===     struct array: contains current point neighbours, see neighbourSearch function header for properties
%   currentPoint    ===     struct: single struct object of point being evaluated, see CBPM function header for properties
%   dx              ===     scalar: space discretisation (the same in x and y directions)
%
%
% Author: Solene Hegarty-Cremer
%%
currentPointInfo = struct('middle', currentPoint.gridPointCoords,...
    'bottomLeft', [currentPoint.gridPointCoords(1) - dx/2, currentPoint.gridPointCoords(2) - dx/2], ...
    'bottomRight', [currentPoint.gridPointCoords(1) + dx/2, currentPoint.gridPointCoords(2) - dx/2], ...
    'topRight', [currentPoint.gridPointCoords(1) + dx/2, currentPoint.gridPointCoords(2) + dx/2], ...
    'topLeft', [currentPoint.gridPointCoords(1) - dx/2, currentPoint.gridPointCoords(2) + dx/2]);

changeOfBasis = [interpPoints(1).normal(2), -interpPoints(1).normal(1);
                  interpPoints(1).normal(1), interpPoints(1).normal(2)];
 
              
localInterpPoints = struct('xCoord', 0, 'yCoord', 0, 'val', interpPoints(1).val, ...
    'vs', interpPoints(1).vs, 'cellType', interpPoints(1).cellType);

for j = 2:length(interpPoints)
   localCoords = changeOfBasis*(interpPoints(j).footPointCoords - interpPoints(1).footPointCoords)';
   localInterpPoints(j) = struct('xCoord', localCoords(1), 'yCoord', localCoords(2), ...
       'val', interpPoints(j).val, 'vs', interpPoints(j).vs, 'cellType', interpPoints(j).cellType);
end

currentPointInfo.middle = changeOfBasis*(currentPointInfo.middle - interpPoints(1).footPointCoords)';
currentPointInfo.bottomLeft = changeOfBasis*(currentPointInfo.bottomLeft - interpPoints(1).footPointCoords)';
currentPointInfo.bottomRight = changeOfBasis*(currentPointInfo.bottomRight - interpPoints(1).footPointCoords)';
currentPointInfo.topLeft = changeOfBasis*(currentPointInfo.topLeft - interpPoints(1).footPointCoords)';
currentPointInfo.topRight = changeOfBasis*(currentPointInfo.topRight - interpPoints(1).footPointCoords)';
end
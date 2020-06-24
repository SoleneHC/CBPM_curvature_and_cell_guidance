function collided = collision(newPoint, newNormal, interpPoints, dx, thetaMin)
%% COLLISION() checks if the given points are imminently colliding
%
% collision uses the normals of the marker particle of interest and its m closest
% points to detect if a collision has occurred or is occurring. This is
% determined by the distance of the marker particles and the angle between
% their associated normal vectors.
%
% INPUTS
%   newPoint        ===     vector(2x1): coordinates of newly resampled marker particle
%   newNormal       ===     vector(2x1): surface normal at newly resampled marker particle
%   interpPoints    ===     struct array: contains current point neighbours, see neighbourSearch function header for properties
%   dx              ===     scalar: space discretisation (the same in x and y directions)
%   thetaMin        ===     scalar: parameter used to detect collisions
%
%
% Author: Solene Hegarty-Cremer
%% 

for i = 1:length(interpPoints)
   if norm(newPoint' - interpPoints(i).footPointCoords,'inf') < dx ...
            && dot(newNormal, interpPoints(i).normal) < cos(thetaMin)
       collided = 1;
       fprintf('Collision Detected, removing point\n')
       return 
   else
       collided = 0;
   end
end
end
function activePoints = CBPM(activePoints, Vx, Vy, fdash, vs, D, dt, dx, x, y, domainLen, ...
    Tend, m, thetaMin, colourAxis, figureNum, stepToPlot, analyticSol, r)
%% CBPM() carries out the cell-based particle method given an initialised surface until Tend
%
% The cell-based particle method (CBPM) is a method for solving PDEs on
% moving interfaces developed by Leung et al. (2009). The CBPM function
% takes in an initialised surface in the correct data structure as well as
% the PDE required to evolve the surface position and the density values.
% There are also two optional inputs for an analytic solution on a circle
% which results in both solutions being plotted on the same axes for
% comparison.
%
% INPUTS
%   activePoints    ===     struct array containing all domain grid cells with properties:
%                               gridPointCoords : coordinates of grid cell center
%                               gridPointIndices : indices of grid cell
%                               footPointCoords : coordinates of marker particle
%                               dist : distance between marker particle and cell center
%                               normal : unit normal vector
%                               curv : local curvature
%                               val : density
%                               newVal : updated density value
%                               newFootPointCoords : updated marker particle coordinates
%                               neighbours : marker particle neighbours
%                               curveLength : length of surface contained in cell
%                               alpha : interpolation coefficients for surface position
%                               beta : interpolation coefficients for density
%                               localPoint : marker particle coordinates in local basis
%                               cellType : type of cell
%                               vs : tangential velocity of cell
%                               zeta : interpolation coefficients for tangential velocity
%                               arcPoint : arc length position
%   Vx              ===     anon. function(point, t): describes change in x position of the interface, takes in point struct and time
%   Vy              ===     anon. function(point, t): describes change in y position of the interface, takes in point struct and time
%   fDash           ===     anon. function(point, t): describes change in density without spatial derivatives, takes in point struct and time
%   vs              ===     anon. function(point, t): describes tangential velocity of cells w.r.t the surface, takes in point struct and time
%   D               ===     scalar: diffusivity of cells
%   dt              ===     scalar: time discretisation
%   dx              ===     scalar: space discretisation (the same in x and y directions)
%   x               ===     vector(1xdomainLen): x locations of grid cell edges
%   y               ===     vector(1xdomainLen): y locations of grid cell edges
%   domainLen       ===     scalar: length of one direction of discretised domain
%   Tend            ===     scalar: time to run simulation to
%   m               ===     vector(1x2): number of neighbours to use for interpolation for surface position and scalar quantities (density) respectively
%   thetaMin        ===     scalar: parameter used to detect collisions
%   colourAxis      ===     vector(1x2): density bounds for colorbar
%
% OPTIONAL INPUTS
%   analyticSol     ===     anon. function(): describes analytic solution for density along a circle of radius r as a function of arc length and time
%   r               ===     anon. function(): describes analytic solution for how the radius changes over time
%
%
% Author: Solene Hegarty-Cremer
%%
%Time loop
t = 1;
%Iterate in time while there are still active cells
while t <= (Tend/dt) && ~isempty(find(~cellfun(@isempty, {activePoints.val}), 1))
    %Movement
    fprintf('Starting Movement\n')
    activePoints = applyMovementAndODE(activePoints, Vx, Vy, fdash, vs, D, t*dt, dt);
    fprintf('Finished Movement\n')
    
    %Resample
    fprintf('Starting Resample\n')
    [activePoints, removalPoints] = interpolateAndResample(activePoints, m, x, y, domainLen, thetaMin, dx);
    fprintf('Finished Resample\n')
    
    %Activation
    fprintf('Starting Activation\n')
    [activePoints, newPoints] = activateNewPoints(activePoints, m, x, y, dx, domainLen, thetaMin, Vx, Vy, t*dt);
    fprintf(['Finished Activation: ' num2str(length(newPoints)) ' points activated\n'])
    
    %Deactivation
    fprintf(['Starting Deactivation: ' num2str(length(removalPoints)) ' points deactivated\n'])
    activePoints = deactivatePoints(activePoints, removalPoints, dx);
    fprintf('Finished Deactivation\n')
    
    
    %Find approximate pore center
    allPoints = [activePoints.footPointCoords];
    poreCenter(1) = (min(allPoints(1:2:end)) + max(allPoints(1:2:end)))/2;
    poreCenter(2) = (min(allPoints(2:2:end)) + max(allPoints(2:2:end)))/2;
    
    
    %Plot at intervals
    if mod(t, stepToPlot) == 0 %89
        figureTitle = ['t = ', num2str(round(t*dt, 2))];
        if nargin < 19
            plotSurface(figureNum, x, y, activePoints, figureTitle, ...
                poreCenter, colourAxis);
            drawnow()
        else
            plotAnalyticComparison(figureNum, x, y, activePoints, figureTitle, ...
                poreCenter, t, analyticSol, dt, r, colourAxis)
            drawnow()
        end
        
    end
    t = t+1;
end

%Plot at final time
figureTitle = ['t = ', num2str((t-1)*dt)];
if nargin < 19
    plotSurface(figureNum, x, y, activePoints, figureTitle, ...
        poreCenter, colourAxis);
    drawnow()
else
    plotAnalyticComparison(figureNum, x, y, activePoints, figureTitle, poreCenter,...
        t-1, analyticSol, dt, r, colourAxis)
    drawnow()
end
hold off






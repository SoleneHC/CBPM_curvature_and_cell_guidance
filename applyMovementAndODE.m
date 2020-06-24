function activePoints = applyMovementAndODE(activePoints, Vx, Vy, fdash, vs, D, t, dt)
%% APPLYMOVEMENTANDODE() applies the specified velocity field to the surface and the ODE to the density values
%
% applyMovementAndODE uses Vx and Vy to evolve the location of the surface
% using Euler's method. fDash is used to evolve the density values using
% forward Euler, fDash should not include any spatial derivatives. The
% transport operator function is called from here.
%
% INPUTS
%   activePoints    ===     struct array: contains all domain grid cells, see CBPM function header for properties
%   Vx              ===     anon. function(point, t): describes change in x position of the interface, takes in point struct and time
%   Vy              ===     anon. function(point, t): describes change in y position of the interface, takes in point struct and time
%   fDash           ===     anon. function(point, t): describes change in density without spatial derivatives, takes in point struct and time
%   vs              ===     anon. function(point, t): describes tangential velocity of cells w.r.t the surface, takes in point struct and time
%   D               ===     scalar: diffusivity of cells
%   t               ===     scalar: current simulation time
%   dt              ===     scalar: time discretisation
%
%
% Author: Solene Hegarty-Cremer
%%

%Get all new values and coords
for i = find(~cellfun(@isempty, {activePoints.val})) %over active points
    activePoints(i).newFootPointCoords = activePoints(i).footPointCoords + ...
       real([Vx(activePoints(i), t), Vy(activePoints(i), t)]).*dt;
    
    % Non surface operation
    activePoints(i).val = activePoints(i).val + fdash(activePoints(i), t)*dt;
    
    %Update info
    activePoints(i).dist = norm(activePoints(i).newFootPointCoords - activePoints(i).gridPointCoords, 'inf');
end

% Transport operator
for i = find(~cellfun(@isempty, {activePoints.val})) %over active points
    % Surface transport operator
    activePoints(i).newVal = applyTransport(activePoints(i), dt, D, t);
end

%Update values
for i = 1:length(activePoints)
    activePoints(i).footPointCoords = real(activePoints(i).newFootPointCoords);
    activePoints(i).val = activePoints(i).newVal;
end

end
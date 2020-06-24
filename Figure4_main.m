%% Figure4_main.m
%
% This script generates the subfigures in Figure 4 of 'INSERT PAPER TITLE'.
% A solution to Equation (eq num) is obtained with the constant simType
% determining which of the two subfigures (no tangential cell velocity,
% constant cell velocity) is solved. The solution is compared with the
% exact solution given by Equation (eq num)
%
%
%
%
%
%
% Author: Solene Hegarty-Cremer

clear activePoints
close all
%% Set simulation constants
%Set which figure to draw
simType = 'b';

%Discretisation and domain
dx = 0.01; %spatial discretisation step
xBounds = [-.6 .6]; yBounds = [-.6 .6]; %spatial domain bounds
x = xBounds(1):dx:xBounds(2); y = yBounds(1):dx:yBounds(2); %discretised spatial domain
domainLen = length(x);
dt = 0.025; %time discretisation step
Tend = 10; %final time
DeltaT = 1; %time interval at which to plot
stepToPlot = DeltaT/dt;

%Initial condition
r0 = 0.25; %initial pore radius
rhoInit = @(s) 0.5*((-3*pi/8<s && s < -pi/8) || (pi/8<s && s<3*pi/8)); %initial density as function of arclength
colourAxis = [-0.2, 0.8]; %min and max values for plotting density

%Equation parameters
c = 0.035; %surface normal velocity
vsMag = 1/10; %magnitude of tangential cell velocity
D = 0; %diffusivity (mm^2/day)

%Parameters for CBPM
m = [3, 3]; %number of local neighbours used for surface and density values, respectively
thetaMin = 5*pi/8; %minimum difference in normal angle to ensure correct neighbours

%Exact Solution
rFunc = @(t) c*t + r0;

%% Functions used to evolve interface and values
% Interface movement
un = @(point) c; 

if simType == 'a' %no tangential velocity
    vs = @(point, t) 0;
    analyticSol = @(s,t) rhoInit(s*exp(t*0))*r0*exp(t*0)/(c*t + r0);
elseif simType == 'b' % constant tangential velocity
    vs = @(point, t) rFunc(t)*(vsMag)*(point.footPointCoords(1) > 0)*((point.footPointCoords(2)>0)*acos((point.footPointCoords(1))/rFunc(t)) ...
        - (point.footPointCoords(2)<0)*acos((point.footPointCoords(1))/rFunc(t)) );
    analyticSol = @(s,t) rhoInit(s*exp(t*vsMag))*r0*exp(t*vsMag)/(c*t + r0);
else
    error('Please choose simType = a or b')
end
   
%Velocity field for interface (Equation SI ???)
Vx = @(point, t) -un(point)*point.normal(1) - vs(point,t)*point.normal(2);
Vy = @(point, t) -un(point)*point.normal(2) + vs(point,t)*point.normal(1);

%ODE for density
fdash = @(point, t) -point.val*un(point)/(un(point)*t + r0) + (simType == 'b')*vsMag*point.val;

%% Initialise surface
%Explicit parameterisation
fx = @(s) r0*cos(s);
fy = @(s) r0*sin(s);
% Derivative of parameterisation for unit normal calculation
dfx = @(s) -r0*sin(s);
dfy = @(s) r0*cos(s);
% Function to be minimised to find closest marker particle to Eulerian grid point
minimisingFunc = @(s, coords) norm([fx(s), fy(s)] - coords, 'inf');

%Initialise empty struct
activePoints(length(x)*length(y)) = struct('gridPointCoords', [], 'gridPointIndices', [], ...
                'footPointCoords', [], 'dist', [], 'normal', [], ...
                'curv', [], 'val', [], 'newVal', [], 'newFootPointCoords', [],...
                'neighbours', [], 'curveLength', [], 'alpha', [], 'beta', [],...
                'localPoint', [], 'cellType', [], 'vs', [], 'zeta', [], 'arcPoint', []);

for i = 1:length(x)
    for j = 1:length(y)
        %Find closest marker particle to grid point
        [sP, dist] = fminbnd(@(s) minimisingFunc(s, [x(i), y(j)]), -pi, pi+pi/2);
        %If marker particle is within the grid cel
        if dist < dx/2
            %Calculate normal
            normal = [-dfy(sP), dfx(sP)];
            normal = normal/norm(normal);
            
            curvature = 1/r0;
            
            value = rhoInit(sP);
          
            %Activate cell
            activePoints(i + domainLen*(j-1)) = struct('gridPointCoords', [x(i), y(j)], 'gridPointIndices', [i, j], ...
                'footPointCoords', [fx(sP), fy(sP)], 'dist', dist, 'normal', normal, ...
                'curv', curvature, 'val', value, 'newVal', value, 'newFootPointCoords', [fx(sP), fy(sP)], ...
                'neighbours', [], 'curveLength', [],'alpha', [], 'beta', [],...
                'localPoint', [], 'zeta', vsMag, 'vs', vsMag*r0*sP, 'arcPoint', r0*sP, 'cellType', 1);
        end
        
    end
end
figureTitle = 'Initial Condition';
[activePoints, removalPoints] = interpolateAndResample(activePoints, m, x, y, domainLen, thetaMin, dx);
plotAnalyticComparison(1, x, y, activePoints, figureTitle, [0,0], 0,...
           analyticSol, dt, rFunc, colourAxis)

%% CBPM
activePoints = CBPM(activePoints, Vx, Vy, fdash, vs, D, dt, dx, x, y, domainLen, ...
    Tend, m, thetaMin, colourAxis, 1, stepToPlot, analyticSol, rFunc);

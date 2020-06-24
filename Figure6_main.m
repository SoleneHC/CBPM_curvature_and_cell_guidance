%% Figure6_main.m
%
% This script generates the subfigures in Figure 5 of 'INSERT PAPER TITLE'.
% A solution to Equation (eq num) is obtained with the constant simType
% determining with of the three subfigures ('a' no tangential cell velocity,
% 'b' constant cell velocity, 'c' changing cell velocity) is solved.
%
%
%
%
%
%
% Author: Solene Hegarty-Cremer

%Clear previous sim data
clear activePoints
close all
%% Set simulation constants
%Set which figure to draw
simType = 'a';

%Discretisation and domain
dx = 0.001; %spatial discretisation step
xBounds = [-0.08 0.08]; yBounds = [-0.08 0.08]; %spatial bounds
x = xBounds(1):dx:xBounds(2); y = yBounds(1):dx:yBounds(2); %discretised domain
domainLen = length(x); 
dt = 0.075; %time discretisation step
Tend = 24; %final time
DeltaT = 3; %time interval at which to plot
stepToPlot = DeltaT/dt;

%Initial condition
r0 = 0.075; %initial pore radius
rhoInit = @(l) 160; %initial density as function of arclength
colourAxis = [64 384]; %min and max values for plotting density

%Equation parameters
k = 7.8125e-06; %secretory rate (cells/mm/day)
vsMag = 0.0025; %magnitude of tangential cell velocity
D = 1e-5; %diffusivity (mm^2/day)

%Parameters for CBPM
m = [3, 5]; %number of local neighbours used for surface and density values, respectively
thetaMin = 5*pi/8; %minimum difference in normal angle to ensure correct neighbours

%% Functions used to evolve interface and values
% Interface movement
un = @(point) -point.val*k; %Equation 1

if simType == 'a' %no tangential velocity
    vs = @(point, t) 0;
elseif simType == 'b' % constant tangential velocity
    vs = @(point, t) vsMag;
elseif simType == 'c' %changing tangential velocity
    vs = @(point, t) vsMag*(t<12.5) - vsMag*(t>=12.5);
else
    error('Please choose simType = a,b, or c')
end
   
%Velocity field for interface (Equation SI ???)
Vx = @(point, t) -un(point)*point.normal(1) - vs(point,t)*point.normal(2);
Vy = @(point, t) -un(point)*point.normal(2) + vs(point,t)*point.normal(1);

%ODE for density
fdash = @(point, t) -(point.val>0)*(point.val * point.curv * un(point));


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
        %If marker particle is within the grid cell
        if dist < dx/2
            %Calculate normal
            normal = [-dfy(sP), dfx(sP)];
            normal = normal/norm(normal);
            
            curvature = 1/r0;
            
            value = rhoInit(sP);
            
            %'Stain' cell types            
            if ((mod(sP,2*pi)) >= pi/8 && (mod(sP,2*pi)) < 2*pi/8)...
                 || ((mod(sP,2*pi)) >= 3*pi/8 && (mod(sP,2*pi)) < 4*pi/8)...
                 || ((mod(sP,2*pi)) >= 5*pi/8 && (mod(sP,2*pi)) < 6*pi/8)...
                 || ((mod(sP,2*pi)) >= 7*pi/8 && (mod(sP,2*pi)) < 8*pi/8)...
                 || ((mod(sP,2*pi)) >= 9*pi/8 && (mod(sP,2*pi)) < 10*pi/8)...
                 || ((mod(sP,2*pi)) >= 11*pi/8 && (mod(sP,2*pi)) < 12*pi/8)...
                 || ((mod(sP,2*pi)) >= 13*pi/8 && (mod(sP,2*pi)) < 14*pi/8)...
                 || ((mod(sP,2*pi)) >= 15*pi/8 && (mod(sP,2*pi)) < 16*pi/8)
                cellType = 1;
            else
                cellType = 0;
                
            end
            
            %Initialise marker particle
            activePoints(i + domainLen*(j-1)) = struct('gridPointCoords', [x(i), y(j)], 'gridPointIndices', [i, j], ...
                'footPointCoords', [fx(sP), fy(sP)], 'dist', dist, 'normal', normal, ...
                'curv', curvature, 'val', value, 'newVal', value, 'newFootPointCoords', [fx(sP), fy(sP)], ...
                'neighbours', [], 'curveLength', [],'alpha', [], 'beta', [],...
                'localPoint', [], 'zeta', [], 'vs', 0, 'cellType', cellType, 'arcPoint', r0*sP);
        end
        
    end
end
figureTitle = 'Initial Condition';
[activePoints, removalPoints] = interpolateAndResample(activePoints, m, x, y, domainLen, thetaMin, dx);
plotSurface([1,2], x, y, activePoints, figureTitle, [0,0], colourAxis);

%% CBPM
activePoints = CBPM(activePoints, Vx, Vy, fdash, vs, D, dt, dx, x, y, domainLen, ...
    Tend, m, thetaMin, colourAxis, [1,2], stepToPlot);

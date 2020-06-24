function [lDomain, densities] = plotAnalyticSolution(analyticSolution, t, r, colourAxis)
%% PLOTANALYTICSOLUTION plots the density over the surface of the analytic solution
%
% plotAnalyticSolution takes in a 1D analytic solution (density over arc
% length) and translates that into 2D for plotting. The discretised arc
% length domain and densities at these points are also output for a 1D
% plot.
%
% INPUTS
%   analyticSol     ===     anon.func(l,t): analytic solution as function of arc length and time
%   t               ===     scalar: current simulation time
%   r               ===     scalar: current pore radius
%   colourAxis      ===     vector(1x2): density bounds for colorbar
%
%
% Author: Solene Hegarty-Cremer
%%
%Get make mapping between l domain and surface position
L = pi;
lDomain = linspace(-L+0.001, L-0.001, 500); 

xCoords = zeros(length(lDomain),1)';
yCoords = zeros(length(lDomain),1)';
densities = zeros(length(lDomain),1)';

for i = 1:length(lDomain)
    xCoord = r*cos(lDomain(i));
    yCoord = r*sin(lDomain(i));
    
    rho = analyticSolution(lDomain(i), t);
    
    xCoords(i) = xCoord;
    yCoords(i) = yCoord;
    densities(i) = rho;    
end

% Plot
hold on
h1 = surface([xCoords;xCoords], [yCoords; yCoords], [zeros(size(xCoords));zeros(size(xCoords))], ...
    [densities; densities], 'EdgeColor', 'interp', 'MarkerFaceColor', 'flat');
set(h1, 'linewidth', 1.5);

xlabel('x')
ylabel('y')
title(['Analytic Solution at t = ' num2str(t)])
colormap(flipud(pink))
caxis(colourAxis)
h = colorbar;
ylabel(h, '\rho')
set(gcf,'color','w');
axis equal

end
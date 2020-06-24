function plotSurface(figureNum, x, y, activePoints, figureTitle, center, colourAxis)
%% PLOTSURFACE() plots the current position of the surface 
%
% plotSurface plots the density along the current surface. If two figure
% numbers are given, the second figure with plot the cell types along the
% surface.
%
% INPUTS
%   figureNum       ===     scalar or vector(1x2): number of figure to plot on
%   x               ===     vector(1xdomainLen): x locations of grid cell edges
%   y               ===     vector(1xdomainLen): y locations of grid cell edges
%   activePoints    ===     struct array: contains all domain grid cells, see CBPM function header for properties
%   figureTitle     ===     string: title to display on figure
%   center          ===     vector(1x2): center of surface
%   colourAxis      ===     vector(1x2): density bounds for colorbar
%
% Author: Solene Hegarty-Cremer
%% 
fprintf('Starting Plot \n')
fig = figure(figureNum(1));
hold on
axis equal 

xPlot = zeros(1, length(find(~cellfun(@isempty, {activePoints.val}))));
yPlot = zeros(size(xPlot));
zPlot = zeros(size(xPlot));
cPlot = zeros(size(xPlot));
angles = zeros(size(xPlot));
j = 1;

if length(figureNum) > 1
    cellTypePlot = zeros(size(xPlot));
end

for i = find(~cellfun(@isempty, {activePoints.val})) %over active points   
    hold on
    
    xPlot(j) = activePoints(i).footPointCoords(1);
    yPlot(j) = activePoints(i).footPointCoords(2);
    cPlot(j) = activePoints(i).val;
    angles(j) = atan2(yPlot(j)-center(2), xPlot(j)-center(1));
    
    if length(figureNum) > 1
        cellTypePlot(j) = round(activePoints(i).cellType);
    end   
    
    j = j + 1;
end

[~, perm] = sort(angles);
xPlot = xPlot(perm);
yPlot = yPlot(perm);
cPlot = cPlot(perm);
xPlot(end+1) = xPlot(1);
yPlot(end+1) = yPlot(1);
cPlot(end+1) = cPlot(1);
zPlot(end+1) = 0;

h1 = surface([xPlot;xPlot], [yPlot; yPlot], [zPlot;zPlot], [cPlot; cPlot], 'EdgeColor', 'interp', 'MarkerFaceColor', 'flat');
set(h1, 'linewidth', 1.5);

colormap(flipud(pink))
caxis(colourAxis)
h = colorbar;
ylabel(h, '\rho')
title(figureTitle)
set(gcf,'color','w');
axis equal
xlim([min(x) max(x)]), ylim([min(y) max(y)])
xlabel('x'), ylabel('y')

if length(figureNum) > 1
    figure(figureNum(2))
    cellTypePlot = cellTypePlot(perm);
    cellTypePlot(end+1) = cellTypePlot(1);
    
    h2  = surface([xPlot;xPlot], [yPlot; yPlot], [zPlot;zPlot], [cellTypePlot; cellTypePlot],...
        'EdgeColor', 'interp', 'MarkerFaceColor', 'flat');
    set(h2, 'linewidth', 1.5);
    
    title(figureTitle)
    set(gcf,'color','w');
    axis equal
    xlim([min(x) max(x)]), ylim([min(y) max(y)])
    xlabel('x'), ylabel('y')
end

drawnow()

fprintf('Plot finished\n')
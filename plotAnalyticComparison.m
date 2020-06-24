function plotAnalyticComparison(figureNum, x, y, activePoints, figureTitle, center, t,...
           analyticSol, dt, rFunc, colourAxis)
%% PLOTANALYTICCOMPARISON() plots the current position of the surface and the analytic solution to the problem
%
% plotAnalyticSolutions plots the density along the current surface. Two
% figures are plotted, one with two subplots with the density along the
% surface and the other with the density along the arc length.
%
% INPUTS
%   figureNum       ===     scalar or vector(1x2): number of figure to plot on
%   x               ===     vector(1xdomainLen): x locations of grid cell edges
%   y               ===     vector(1xdomainLen): y locations of grid cell edges
%   activePoints    ===     struct array: contains all domain grid cells, see CBPM function header for properties
%   figureTitle     ===     string: title to display on figure
%   center          ===     vector(1x2): center of surface
%   t               ===     scalar: current simulation timestep number
%   analyticSol     ===     anon.func(l,t): analytic solution as function of arc length and time
%   dt              ===     scalar: timestep length
%   rFunc           ===     anon.func(l,t): analytic solution for pore radius as function of time
%   colourAxis      ===     vector(1x2): density bounds for colorbar
%
% Author: Solene Hegarty-Cremer          
%% 
r = rFunc(t*dt);

%Plot whole surface, two subplots - one for analytic, one for CBPM
fig1 = figure(figureNum(1));
hold on
subplot(1,2,2)
[paraX, uX] = plotAnalyticSolution(analyticSol, t*dt, r, colourAxis);
xlim([min(x) max(x)])
ylim([min(y) max(y)])

subplot(1,2,1)
plotSurface(figureNum, x, y, activePoints, figureTitle, center, colourAxis);

%Solution over arc length
fig2 = figure(max(figureNum)+1);
colours = lines; 
cl = colours(round(t*dt) + 1,:);

hold on
plot(paraX, uX, 'Color', cl, 'LineWidth', 1, 'DisplayName', ['t = ' num2str(round(t*dt,1))]);
for i = find(~cellfun(@isempty, {activePoints.val}))
    r = sqrt(activePoints(i).footPointCoords(1)^2 + activePoints(i).footPointCoords(2)^2);
    if activePoints(i).footPointCoords(2) > 0
        s = acos((activePoints(i).footPointCoords(1))/r);
    else
        s = -acos((activePoints(i).footPointCoords(1))/r);
    end

    scatter(s, activePoints(i).val, 10, cl, 'LineWidth', 1, 'HandleVisibility','off')
    hold on
end
xlabel('l')
ylabel('\rho')
ylim(colourAxis)
xlim([-pi pi])
legend show

drawnow()

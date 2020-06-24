function [newLocalPoint, curvature, newNormalLocal, dist, neighbours, ...
    curveLength, toRemove] = localResample(alpha, localGridPoint)
%% LOCALRESAMPLE Resamples surface given local interpolation
%
% localResample uses the current grid cell rotated into local coordinate
% and checks whether the interpolated surface passes through the cell. If
% yes, the point on the local interpolation closest to the cell center is
% chosen as the new marker particle location in local coordinates
%
% INPUTS
%   alpha               ===     vector(3x1): describes parabola fit through local points
%   localGridPoint      ===     struct: describes current grid cell in local coordinates, with properties
%                                   middle
%                                   bottomLeft
%                                   bottomRight
%                                   topLeft
%                                   topRight
%
%
% Author: Solene Hegarty-Cremer
%%
%Initialise local interpolation function
fLocal = @(x) alpha(1) + alpha(2)*x + alpha(3)*x^2;

neighbours = [];
intersections = [];

%Pull edge names from local grid point
edges = fieldnames(localGridPoint);

%Loop through edges to check if fLocal intersects with them
for i = 2:length(edges)
    %Middle of array
    if i <= length(edges)-1
        vert1 = edges{i};
        vert2 = edges{i+1};
    %End of array loop through to start
    else
        vert1 = edges{i};
        vert2 = edges{2};
    end
    
    %Pull coordinates from vertices
    xs = [localGridPoint.(vert1)(1) localGridPoint.(vert2)(1)];
    ys = [localGridPoint.(vert1)(2) localGridPoint.(vert2)(2)];
    
    %If xs form a vertical line
    if abs(localGridPoint.(vert2)(1) - localGridPoint.(vert1)(1)) < 1e-5 %vertical line
        %If flocal passes through vertical line
        if fLocal(localGridPoint.(vert2)(1)) > min(ys) && fLocal(localGridPoint.(vert2)(1)) < max(ys)
            xRoots = localGridPoint.(vert2)(1);
        else %no xRoots
            xRoots = 1i;
        end
    else
        %Make a line between vertices (y=mx_c)
        m = (localGridPoint.(vert2)(2) - localGridPoint.(vert1)(2))/(localGridPoint.(vert2)(1) - localGridPoint.(vert1)(1));
        c = localGridPoint.(vert2)(2) - m*localGridPoint.(vert2)(1);
        % If surface interpolation is a parabola
        if alpha(3) > 1e-10
            xRoots = [(-(alpha(2)-m) + sqrt((alpha(2)-m)^2 - 4*alpha(3)*(alpha(1)-c)))/(2*alpha(3)), ...
                (-(alpha(2)-m) - sqrt((alpha(2)-m)^2 - 4*alpha(3)*(alpha(1)-c)))/(2*alpha(3))];
        else %if surface interpolation is a straight line
            xRoots = (alpha(1)-c)/(m - alpha(2));
        end
    end
    
    %If intersections are real
    for j = 1:length(xRoots)
        if isreal(xRoots(j)) && xRoots(j) >= min(xs) && xRoots(j) <= max(xs)...
                && (abs(xRoots(j)) > 1e-10)
            intersections = [intersections xRoots(j)];
            
            %Get neighbours
            if contains(vert1, 'Left') && contains(vert2, 'Left')
                neighbours = [neighbours; -1 0];
            elseif contains(vert1, 'Right') && contains(vert2, 'Right')
                neighbours = [neighbours; 1 0];
            elseif contains(vert1, 'bottom') && contains(vert2, 'bottom')
                neighbours = [neighbours; 0 -1];
            elseif contains(vert1, 'top') && contains(vert2, 'top')
                neighbours = [neighbours; 0 1];
            end
        end
    end
end

if length(intersections) < 2 %if flocal does not intersect with grid cell
    toRemove = 1;
    dist = 100;
    newFootPointX = 0;
    curveLength = 0;
else %find closest point to cell center
    statPoint = (-alpha(2)/(2*alpha(3)));
    newXoptions = sum(intersections)/length(intersections);
    
    if length(intersections) > 2
        fprintf('More than two intersections \n')
        sorted = sort(intersections);
        intersections = sorted(1:2);
        newXoptions = sum(intersections)/length(intersections);
        
    elseif statPoint >= min(xs) && statPoint <= max(xs) && fLocal(statPoint) > min(ys) && fLocal(statPoint) < max(ys)
        newXoptions = [newXoptions statPoint];
    end
    
    %Find distance between cell center and x options
    infNorms = zeros(length(newXoptions),1);
    for i = 1:length(infNorms)
        infNorms(i) = norm(localGridPoint.middle' - [newXoptions(i), fLocal(newXoptions(i))], 'inf');
    end
    
    [dist, bestOption] = min(infNorms);
    newFootPointX = newXoptions(bestOption);
    toRemove = 0;

    %Estimate length of curve contained in grid cell
    curveLength = abs(intersections(2) - intersections(1));
end

%Get coords
newLocalPoint = [newFootPointX; fLocal(newFootPointX)];

%Get new normal
oldNormalLocal = [0; 1];
newNormalLocal = getNewNormal(newLocalPoint, localGridPoint, oldNormalLocal, alpha);

%Update geometry scalars
curvature = 2*alpha(3)/(1 + (alpha(2) + 2*alpha(3)*newFootPointX)^2)^(3/2);
end
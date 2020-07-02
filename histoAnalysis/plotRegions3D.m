function plotRegions3D(xyzCoords, size, color, xrange, yrange, transparency)

% This function is to plot the 3-dimension logical matrix (created by
% reformatTiffStack function and getCoordinates function) in 3-d space.

% INPUT: 



if ~exist('size', 'var'); size = 10; end
%if ~exist('c', 'var'); color = 'b'; end
if ~exist('xrange', 'var'); xrange = 4000; end
if ~exist('yrange', 'var'); yrange = 4000; end
if ~exist('transparency', 'var'); transparency = 1; end

h = scatter3(xyzCoords(:, 1), xyzCoords(:, 2), xyzCoords(:, 3), size, color, 'filled');
set(h, 'MarkerEdgeAlpha', transparency, 'MarkerFaceAlpha', transparency);

xlim([xrange(1), xrange(2)]);
ylim([yrange(1), yrange(2)]);
set( gca, 'ZDir', 'reverse' );
set(gca, 'XDir', 'normal');
set(gca, 'YDir', 'normal');
xlabel('ML axis (in micron)')
ylabel('AP axis (in micron)')
zlabel('DV axis (in micron)')

end
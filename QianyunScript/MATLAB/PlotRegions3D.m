function PlotRegions3D(xyzCoords, size, color, xrange, yrange)

% This function is to plot the 3-dimension logical matrix (created by
% reformatTiffStack function and getCoordinates function) in 3-d space.

% INPUT: 



if ~exist('size', 'var'); size = 10; end
%if ~exist('c', 'var'); color = 'b'; end
if ~exist('xrange', 'var'); xrange = 2000; end
if ~exist('yrange', 'var'); yrange = 1000; end

scatter3(xyzCoords(:, 1), xyzCoords(:, 3), xyzCoords(:, 2), size, color, 'filled');
xlim([0, xrange]);
ylim([0, yrange]);
set( gca, 'ZDir', 'reverse' );
set(gca, 'XDir', 'normal');
set(gca, 'YDir', 'normal');
xlabel('ML axis (in pixel)')
ylabel('AP axis (in pixel)')
zlabel('DV axis (in pixel)')

end
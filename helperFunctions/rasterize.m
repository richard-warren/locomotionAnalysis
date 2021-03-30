function rasterize(ax)
% rasterize an axis while keeping axis labels intact // doesn't support 3D plots

if nargin==0; ax = gca; end
xlims = xlim;
ylims = ylim;


% get image of current axis (just the axis contents - not the labels and edges, etc.)
set(gca, 'visible', 'off')  % make sure axis edges aren't included in image of axis
img = getframe(ax);
img = flipud(frame2im(img));
set(gca, 'visible', 'on')

cla  % remove existing objects that would otherwise be vectorized
image(linspace(xlims(1), xlims(2), size(img,2)), ...
    linspace(ylims(1), ylims(2), size(img,1)), ...
    img)

% make sure no funny business changes to limits, and make sure axis remains on top of timage
set(gca, 'xlim', xlims, 'ylim', ylims, 'layer', 'top')
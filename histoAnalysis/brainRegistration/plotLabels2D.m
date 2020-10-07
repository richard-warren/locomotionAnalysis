function plotLabels2D(labels, varargin)


% settings
s.dim = 'ap';  % which dimension to collapse across
s.apGrid = 1:size(labels,1);
s.mlGrid = 1:size(labels,3);
s.dvGrid = 1:size(labels,2);
s.colors = [];
s.patchArgs = {};
s.smoothing = 1;                  % (samples) smoothing for traces


% inits
regionLabels = unique(labels(:));
regionLabels = regionLabels(regionLabels~=0)';
s.colors = lines(length(unique(labels(:)))-1);
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
hold on; axis equal


% view specific settings
switch s.dim
    case 'ap'
        projectionDim = 1;
        axes = {s.mlGrid, s.dvGrid};
        xdir = 'normal';
        ydir = 'reverse';
        axLabels = {'ML', 'DV'};
        transpose = true;
    
    case 'dv'
        projectionDim = 2;
        axes = {s.mlGrid, s.apGrid};
        xdir = 'normal';
        ydir = 'normal';
        axLabels = {'ML', 'AP'};
        transpose = true;
        
    case 'ml'
        projectionDim = 3;
        axes = {s.apGrid, s.dvGrid};
        xdir = 'reverse';
        ydir = 'reverse';
        axLabels = {'AP', 'DV'};
        transpose = false;
end



for i = regionLabels
    projection = squeeze(max(labels==i, [], projectionDim));
    if transpose; projection = projection'; end
%     projection = squeeze(mean(labels==i, projectionDim));
%     projection = projection > (max(projection(:)) * .5);
    
    [i1, i2] = ind2sub(size(projection), find(projection));
    ax1 = axes{1}(i1);
    ax2 = axes{2}(i2);
    k = boundary(ax1', ax2');
    x = ax1(k)';
    y = ax2(k)';
    
    % interpolate via arc-length parameterization
    x = smooth(x, s.smoothing);
    y = smooth(y, s.smoothing);
    sl = [0 cumsum(sqrt(sum(diff([x'; y'],1,2).^2,1)))];  % segment length
    
    % interpolation parameterized by arc length :)
    method = 'maxima';
    x = interp1(sl, x, linspace(sl(1), sl(end), 200), method);
    y = interp1(sl, y, linspace(sl(1), sl(end), 200), method);
    
    
    patch(x, y, s.colors(i,:), 'FaceAlpha', .4, 'LineWidth', 2, 'EdgeColor', s.colors(i,:), s.patchArgs{:});
end

set(gca, 'ydir', ydir, 'xdir', xdir);
xlabel(axLabels{1})
ylabel(axLabels{2})

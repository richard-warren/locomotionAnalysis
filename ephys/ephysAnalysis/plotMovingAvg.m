function plotMovingAvg(x, y, varargin)

% given two continuous variables, x and y, plots moving average of y as a
% function of x


% settings
% --------

% moving average
s.percentileLims = [.1 99.9];  % (percentiles) x axis limits
s.binNum = 500;  % number of points on the x axis
s.windowSize = .05;  % (fraction of x axis) width of moving average window
s.smoothing = .05;  % (fraction of x axis) width of moving average smoothing

% plot
s.scatters = 1000;  % max number of scatter points to show
s.color = [.15 .15 .15];
s.ylabel = '';
s.xlabel = '';
s.plotDensity = true;  % whether to plot pdf in background
s.colorDensity = [.4 .4 1];;


% initializations
% ---------------
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end % reassign settings passed in varargin
x = x(:); y = y(:);  % make columns
xLims = prctile(x, s.percentileLims);
xGrid = linspace(xLims(1), xLims(2), s.binNum);
winSz = s.windowSize * diff(xLims);
yBins = ~isnan(y);
s.scatters = min(s.scatters, sum(yBins));

[movMean, movStd, density] = deal(nan(1, length(xGrid)));
for i = 1:length(xGrid)
    xBins = x>=(xGrid(i)-.5*winSz) & x<(xGrid(i)+.5*winSz);
    density(i) = sum(xBins);
    
    bins = xBins & yBins;
    if any(bins)
        movMean(i) = sum(y(bins)) / sum(bins);
        movStd(i) = std(y(bins));
    end
end
density = density / sum(density);  % normalize


% smooth
if s.smoothing
    smoothSmps = round((s.smoothing * range(xLims)) / (xGrid(2)-xGrid(1)));
    smoothSmps = smoothSmps + double(mod(smoothSmps,2)==0);  % ensure odd number
    
    movMean = smooth(movMean, smoothSmps)';
    movStd = smooth(movStd, smoothSmps)';
    density = smooth(density, smoothSmps)';
end


% plot
% ----

close all
figure('color', 'white', 'menubar', 'none'); hold on

% error bars
patch([xGrid(1) xGrid fliplr(xGrid)], ...
    [movMean(1)-movStd(1) movMean+movStd fliplr(movMean-movStd)], ...
    s.color, 'EdgeColor', 'none', 'FaceAlpha', .2)

% mean
plot(xGrid, movMean, 'Color', s.color, 'LineWidth', 2)

% scatter
xy = [x y];
xy = datasample(rmmissing(xy), s.scatters);
scatter(xy(:,1), xy(:,2), 10, s.color, 'filled', ...
    'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .4)
yLims = get(gca, 'YLim');

% density
if s.plotDensity
    densityScaled = density * (range(yLims) / max(density)) - -yLims(1);
    patch([xGrid(1) xGrid xGrid(end)], ...
        [yLims(1) densityScaled yLims(1)], ...
        s.colorDensity, 'FaceAlpha', .2, 'EdgeColor', 'none') 
end

% fancify
set(gca, 'xlim', xLims, 'ylim', yLims)
if ~isempty(s.xlabel); xlabel(s.xlabel); end
if ~isempty(s.ylabel); ylabel(s.ylabel); end


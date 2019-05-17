function heatmapRick(x, y, opts)

% plots joint pdf for x and y and plots as heat map

% settings
s.xlabel = [];
s.ylabel = [];
s.xLims = [];
s.yLims = [];
s.percentileLims = [5 95]; % exclude x values outside these percentile limits
s.binNum = 100;

% reassign settings contained in opts
if exist('opts', 'var'); for i = 1:2:length(opts); s.(opts{i}) = opts{i+1}; end; end



% initializations
if isempty(s.xLims); s.xLims = prctile(x, s.percentileLims); end % 5th to 95th percentile
if isempty(s.yLims); s.yLims = prctile(y, s.percentileLims); end % 5th to 95th percentile
xBins = linspace(s.xLims(1), s.xLims(2), s.binNum);
yBins = linspace(s.yLims(1), s.yLims(2), s.binNum);
[meshX, meshY] = meshgrid(xBins, yBins);

% compute heatmap
heatmap =  ksdensity([x(:), y(:)], [meshX(:), meshY(:)]);
heatmap = reshape(heatmap, length(xBins), length(yBins));
heatmap = heatmap ./ repmat(max(heatmap, [], 1), length(yBins), 1); % set max of each column to 1

colormap hot
imagesc(xBins, yBins, heatmap)
set(gca, 'YDir', 'normal', 'TickDir', 'out', 'Box', 'off')

if ~isempty(s.xlabel); xlabel(s.xlabel); end
if ~isempty(s.ylabel); ylabel(s.ylabel); end
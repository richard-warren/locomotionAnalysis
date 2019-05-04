function heatmapRick(x, y, varNames, xLims, yLims)


% settings
percentileLims = [5 95]; % exclude x values outside these percentile limits
bins = 100;

% initializations
if ~exist('xLims', 'var') xLims = prctile(x, percentileLims); end % 5th to 95th percentile
if ~exist('yLims', 'var') yLims = prctile(y, percentileLims); end % 5th to 95th percentile
xBins = linspace(xLims(1), xLims(2), bins);
yBins = linspace(yLims(1), yLims(2), bins);
[meshX, meshY] = meshgrid(xBins, yBins);

% compute heatmap
heatmap =  ksdensity([x(:), y(:)], [meshX(:), meshY(:)]);
heatmap = reshape(heatmap, length(xBins), length(yBins));
heatmap = heatmap ./ repmat(max(heatmap, [], 1), length(yBins), 1); % set max of each column to 1

colormap hot
imagesc(xBins, yBins, heatmap)
set(gca, 'YDir', 'normal', 'DataAspectRatio', [1 1 1], 'TickDir', 'out', 'Box', 'off')
xlabel(varNames{1})
ylabel(varNames{2})

%% SUCCESS BY HEIGHT FOR EACH CONDITION

% settings
bins = 
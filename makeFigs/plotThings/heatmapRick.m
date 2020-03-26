function [heatmap, xBins, yBins] = heatmapRick(x, y, varargin)

% plots joint pdf for x and y and plots as heat map

% settings
s.xlabel = [];
s.ylabel = [];
s.xLims = [];
s.yLims = [];
s.percentileLims = [5 95]; % exclude x values outside these percentile limits
s.binNum = 100;
s.colormap = 'hot';
s.showPlot = true;  % whether to plot heatmap (otherwise you can just use the returned heatmap)
s.normalize = '';  % 'row' or 'col' normalize the sum of each row or column to 1


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings contained in opts
if isempty(s.xLims); s.xLims = prctile(x, s.percentileLims); end % 5th to 95th percentile
if isempty(s.yLims); s.yLims = prctile(y, s.percentileLims); end % 5th to 95th percentile
xBins = linspace(s.xLims(1), s.xLims(2), s.binNum);
yBins = linspace(s.yLims(1), s.yLims(2), s.binNum);
[meshX, meshY] = meshgrid(xBins, yBins);

% compute heatmap
heatmap =  ksdensity([x(:), y(:)], [meshX(:), meshY(:)]);
heatmap = reshape(heatmap, length(xBins), length(yBins));
if strcmp(s.normalize, 'row')
    heatmap = heatmap ./ repmat(max(heatmap, [], 2), 1, length(xBins)); % set max of each row to 1
elseif strcmp(s.normalize, 'col')
    heatmap = heatmap ./ repmat(max(heatmap, [], 1), length(yBins), 1); % set max of each column to 1
end
    

if s.showPlot
    colormap(s.colormap)
    imagesc(xBins, yBins, heatmap)
    set(gca, 'YDir', 'normal', 'TickDir', 'out', 'Box', 'off')
end

if ~isempty(s.xlabel); xlabel(s.xlabel); end
if ~isempty(s.ylabel); ylabel(s.ylabel); end
function [corrs, slopes] = scatterPlotRick(x, y, conditions, opts)

% scatters x vs y and fits a line through them // conditions contains
% condition # for each xy pair // dft conditions are plotted in dft colors
% and each have their own line fit! // returns corrs, the correlation btwn
% x and y for each condition

% settings
s.colors = 'hsv'; % color scheme // can be specified either as a string, or as an nX3 color map, where n is number of conditions
s.scatSize = 20; % scatter size
s.scatAlpha = .2; % scatter circle transparency
s.lineAlpha = 1; % line transparency
s.maxScatterPoints = 2000; % don't show any more than this many scatter points
s.conditionNames = {}; % cell array with names of conditions // user can specify this by passing in via 'opts' // if specified, a legend is added


% reassign settings contained in opts
if exist('opts', 'var'); for i = 1:2:length(opts); s.(opts{i}) = opts{i+1}; end; end

% initializations
x = x(:); y=y(:); conditions=conditions(:); % ensure similar orientation
conditionNum = length(unique(rmmissing(conditions)));
if ischar(s.colors); s.colors = eval([s.colors '(conditionNum)']); end % set colorspace if color is specified as a string
scatters = nan(1, conditionNum);
lines = nan(1, conditionNum);
scatterPointsPerCondition = round(s.maxScatterPoints / conditionNum);
corrs = nan(1,conditionNum);
slopes = nan(1,conditionNum);


% plot everything
for h = 1:conditionNum
    
    % get paw and obs height for bin
    bins = conditions==h & ~isnan(x) & ~isnan(y);
    binsSub = bins;
    if sum(bins)>scatterPointsPerCondition
        inds = find(bins);
        inds = inds(randperm(length(inds), scatterPointsPerCondition));
        binsSub = false(size(conditions));
        binsSub(inds) = true;
    end
    
    % scatter that shit
    scatters(h) = scatter(x(binsSub), y(binsSub), s.scatSize, s.colors(h,:), 'filled', ...
        'MarkerEdgeAlpha', s.scatAlpha, 'MarkerFaceAlpha', s.scatAlpha); hold on
    
    % add best fit line
    fit = polyfit(x(bins), y(bins), 1);
    lines(h) = plot(x(bins), polyval(fit, x(bins)), 'linewidth', 4, 'color', [s.colors(h,:) s.lineAlpha]);
    corrs(h) = corr(x(bins), y(bins));
    slopes(h) = fit(1);
end

% pimp fig
uistack(lines, 'top')
if ~isempty(s.conditionNames); legend(lines, s.conditionNames, 'location', 'best', 'box', 'off'); end




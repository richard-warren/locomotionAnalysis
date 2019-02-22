function logPlotRick(x, y, varNames, conditions, conditionNames)

% plots probability of logical variable being true (y axis) as function of
% scalar (x axis) // x and y are the scalar and logical variable,
% respectively // varNames is a cell array containing names of x, then y
% value // conditions is of length x, indicating the condition identity of
% values is x and y // each condition is plotted as a line of a different
% color // the names of the conditions are stored in conditionNames, a cell
% array of length max(conditions)

% settings
xPercentileLims = [5 95]; % exclude x values outside these percentile limits
smoothing = .1; % expressed as fraction of x limits (determined by xPercentileLims)
bins = 100;
lineWidth = 3;



% initializations
colors = hsv(max(conditions));

xLims = prctile(x, xPercentileLims); % 5th to 95th percentile
sigma = range(xLims) * smoothing; % .1 of standard deviation
binEdges = min(x) : diff(xLims)/bins : max(x);
binCenters = binEdges(1:end-1) + .5*mode(diff(binEdges));

kernel = arrayfun(@(x) (1/(sigma*sqrt(2*pi))) * exp(-.5*(x/sigma)^2), ...
    -sigma*5:mode(diff(binEdges)):sigma*5);
kernel = kernel / sum(kernel);



for i = 1:max(conditions)
    
    binAvgs = histcounts(x(conditions==i & y), binEdges) ./ ...
                histcounts(x(conditions==i), binEdges);
    binAvgs = fillmissing(binAvgs, 'linear');
    histoConv = conv(binAvgs, kernel, 'same');
    plot(binCenters, histoConv, 'LineWidth', lineWidth, 'Color', colors(i,:)); hold on
    
    if exist('subIdentities', 'var')
        for j = unique(subIdentities)'
            binAvgs = histcounts(x(conditions==i & subIdentities==j & y), binEdges) ./ ...
                histcounts(x(conditions==i & subIdentities==j), binEdges);
            binAvgs = fillmissing(binAvgs, 'linear');
            histoConv = conv(binAvgs, kernel, 'same');
            plot(binCenters, histoConv, 'LineWidth', thinLineWidth, 'Color', colors(i,:)); hold on
        end
    end
    
    
    
end

set(gca, 'YLim', [0 1], 'XLim', xLims, 'Box', 'off')
xlabel(varNames{1})
ylabel(varNames{2})
if exist('conditionNames', 'var'); legend(conditionNames, 'Box', 'off', 'Location', 'best'); end


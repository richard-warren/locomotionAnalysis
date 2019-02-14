function logPlotRick(x, y, varNames, conditions, conditionNames)

% temp
% flat = getNestedStructFields(data, {'mouse', 'session', 'trial', 'isLightOn', 'obsHgt', ...
%     'isWheelBreak', 'isBigStep', 'condition', 'isTrialSuccess', 'trialVel', 'trialAngle'});
% conditionNames = unique({flat.condition});
% x = [flat.trialVel];
% y = [flat.isLightOn];
% conditions = cellfun(@(x) find(ismember(conditionNames,x)), {flat.condition});
% varNames = {'x', 'y'};
% close all; figure('Position', [1930 552 761 355], 'color', 'white', 'MenuBar', 'none'); 


% settings
xPercentileLims = [5 95]; % exclude x values outside these percentile limits
smoothing = .04; % expressed as fraction of x limits (determined by xPercentileLims)
bins = 100;



% initializations
colors = hsv(length(conditionNames));

xLims = prctile(x, xPercentileLims); % 5th to 95th percentile
sigma = range(xLims) * smoothing; % .1 of standard deviation
binEdges = linspace(xLims(1), xLims(2), bins+1);
binCenters = binEdges(1:end-1) + .5*mode(diff(binEdges));

kernel = arrayfun(@(x) (1/(sigma*sqrt(2*pi))) * exp(-.5*(x/sigma)^2), ...
    -sigma*5:mode(diff(binEdges)):sigma*5);
kernel = kernel / sum(kernel);



for i = 1:length(conditionNames)
    
    binAvgs = histcounts(x(conditions==i & y), binEdges) ./ ...
                histcounts(x(conditions==i), binEdges);
    binAvgs = fillmissing(binAvgs, 'linear');
    histoConv = conv(binAvgs, kernel, 'same');
    
    plot(binCenters, histoConv, 'LineWidth', 2, 'Color', colors(i,:)); hold on
    
end

set(gca, 'YLim', [0 1], 'XLim', xLims, 'Box', 'off')
xlabel(varNames{1})
ylabel(varNames{2})
legend(conditionNames, 'Box', 'off')


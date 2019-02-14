function scatterPlotRick(vars, varNames, conditions, conditionNames)

% settings
colors = hsv(length(conditionNames));
circSize = 20;
transparency = .2;
maxScatterPoints = 2000;


% initializations
scatters = nan(1, length(conditionNames));
lines = nan(1, length(conditionNames));


% plot everything
for h = 1:length(conditionNames)
    
    % get paw and obs height for bin
    conditionBins = conditions==h & all(~isnan(vars),1);
    conditionBinsSub = conditionBins;
    if sum(conditionBins)>maxScatterPoints
        inds = find(conditionBins);
        inds = inds(randperm(length(inds), maxScatterPoints));
        conditionBinsSub = false(size(conditions));
        conditionBinsSub(inds) = true;
    end
        
    
    % scatter that shit
    scatters(h) = scatter(vars(1,conditionBinsSub), vars(2,conditionBinsSub), circSize, colors(h,:), 'filled', ...
        'MarkerEdgeAlpha', transparency, 'MarkerFaceAlpha', transparency); hold on
    
    % add best fit line
    fit = polyfit(vars(1,conditionBins), vars(2,conditionBins), 1);
    lines(h) = plot(vars(1,conditionBins), polyval(fit, vars(1,conditionBins)), 'linewidth', 4, 'color', colors(h,:));
end

% add obs height line
plot(vars(1,conditionBins), vars(1,conditionBins), '-', 'linewidth', 2, 'color', [.5 .5 .5])

% pimp fig
xlabel(varNames{1})
ylabel(varNames{2})
uistack(lines, 'top')
if length(conditionNames)>1; legend(lines, conditionNames, 'location', 'northwest', 'box', 'off'); end
pbaspect([1 1 1])

% blackenFig




function scatterObsVsPawHeights(data, bins, binLabels)


% settings
colors = hsv(length(binLabels));
validBins = ~[data.isWheelBreak];
circSize = 20;
transparency = .2;
xLims = [3.175 10];

% initializations
bins(~validBins) = 0;
close all; figure('position', [1000 250 475 425], 'color', 'white', 'menubar', 'none');

scatters = nan(1, length(binLabels));
lines = nan(1, length(binLabels));

for h = 1:length(binLabels)
    
    % get paw and obs height for bin
    obsHgts = [data(bins==h).obsHeightsVid];
    pawHgts = cell2mat(cellfun(@(x,ind) max(squeeze(x{data(ind).firstPawOver}(end,3,:))), ...
        {data(bins==h).modifiedLocations}, num2cell(find(bins==h)), 'UniformOutput', false)) * 1000;
    
    invalidBins = pawHgts<obsHgts;
    pawHgts(invalidBins) = nan;
    obsHgts(invalidBins) = nan;
    
    % scatter that shit
    scatters(h) = scatter(obsHgts, pawHgts, circSize, colors(h,:), 'filled', ...
        'MarkerEdgeAlpha', transparency, 'MarkerFaceAlpha', transparency); hold on
    
    % add best fit line
    fit = polyfit(obsHgts(~invalidBins), pawHgts(~invalidBins), 1);
    lines(h) = plot(obsHgts, polyval(fit, obsHgts), 'linewidth', 5, 'color', colors(h,:)*.5);
end

% add obs height line
plot(obsHgts, obsHgts, '-', 'linewidth', 2, 'color', [.5 .5 .5])

% pimp fig
xlabel('obstacle height (mm)')
ylabel('paw height (mm)')
uistack(lines, 'top')
if length(binLabels)>1; legend(scatters, binLabels, 'location', 'northwest'); end
set(gca, 'XLim', xLims)
pbaspect([1 1 1])





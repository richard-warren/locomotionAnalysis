function scatterObsVsPawHeights(data, bins, binLabels)


% settings
colors = hsv(length(binLabels))+.2; colors(colors>1)=1;
% validBins = ~[data.isWheelBreak];
validBins = ones(1,length(data));
circSize = 20;
transparency = .2;
xLims = [3.175 10];
preObsLim = .008;

% initializations
figure('position', [1000 250 475 425], 'color', 'white', 'menubar', 'none', 'InvertHardCopy', 'off');

scatters = nan(1, length(binLabels));
lines = nan(1, length(binLabels));


% get all paw heights at the moment paw is preObsLim in front of obs
pawHgts = nan(1,length(data));
for i = 1:length(data)
    firstPawOver = data(i).firstPawOver;
    locations = squeeze(data(i).modifiedLocations{firstPawOver}(end,:,:));
    hgtInd = find(locations(1,:)>-preObsLim,1,'first');
    if hgtInd>1; pawHgts(i) = locations(3, hgtInd)*1000; end
end
obsHgts = [data.obsHeightsVid];

% remove bins where paw is never higher than obs
pawMaxHgts = cell2mat(cellfun(@(x,ind) max(squeeze(x{data(ind).firstPawOver}(end,3,:))), ...
             {data.modifiedLocations}, num2cell(1:length(data)), 'UniformOutput', false)) * 1000;
bins(pawMaxHgts<obsHgts | isnan(obsHgts) | isnan(pawHgts) | ~validBins) = 0;


% plot everything
for h = 1:length(binLabels)
    
    % get paw and obs height for bin
    obsHgtsBin = obsHgts(bins==h);
    pawHgtsBin = pawHgts(bins==h);
    
    % scatter that shit
    scatters(h) = scatter(obsHgtsBin, pawHgtsBin, circSize, colors(h,:), 'filled', ...
        'MarkerEdgeAlpha', transparency, 'MarkerFaceAlpha', transparency); hold on
    
    % add best fit line
    fit = polyfit(obsHgtsBin, pawHgtsBin, 1);
    lines(h) = plot(obsHgtsBin, polyval(fit, obsHgtsBin), 'linewidth', 4, 'color', colors(h,:));
%     lines(h) = plot(obsHgtsBin, polyval(fit, obsHgtsBin), 'linewidth', 4, 'color', 'black');
end

% add obs height line
plot(obsHgtsBin, obsHgtsBin, '-', 'linewidth', 2, 'color', [.5 .5 .5])

% pimp fig
xlabel('obstacle height (mm)')
ylabel('paw height (mm)')
uistack(lines, 'top')
if length(binLabels)>1; legend(lines, binLabels, 'location', 'northwest', 'box', 'off'); end
set(gca, 'XLim', xLims)
pbaspect([1 1 1])

blackenFig




function scatterHistoRick(x, y, opts)

% scatters data with nice marginal histos on the top and right of the plot
% // can give grouped data such that dft groups of data have dft colors,
% along with corresponding dft histos on the sides

% settings
s.colors = 'hsv'; % color scheme // can be specified either as a string, or as an nX3 color map, where n is number of data points
s.xLims = [];
s.yLims = [];
s.scatAlpha = .08;
s.scatPlotSize = .65; % fraction of plot to be occupied by the scatter
s.border = .15; % fraction of plot to be occupied by the border
s.xlabel = [];
s.ylabel = [];
s.groupId = []; % group id for each element in x and y, categorical variable deterimining color of scatter point 
s.groupFcn = @nanmean; % summary state used to collapse across items within a group (eg. change to median if you want to take the median within each group)
s.showCrossHairs = true; % whether to show crosshairs showing std and mean for each variable

% reassign settings contained in opts
if exist('opts', 'var'); for i = 1:2:length(opts); s.(opts{i}) = opts{i+1}; end; end



% initializations
if isempty(s.groupId); s.groupId = ones(size(x)); end
if isempty(s.xLims); s.xLims = [min(x) max(x)]; end
if isempty(s.yLims); s.yLims = [min(y) max(y)]; end
xGrid = linspace(s.xLims(1), s.xLims(2), 200);
yGrid = linspace(s.yLims(1), s.yLims(2), 200);
if ischar(s.colors); s.colors = eval([s.colors '(max(s.groupId))']); end % set colorspace if color is specified as a string
groupNum  = max(s.groupId);

% get summary stats within each group
if groupNum>1
    [groupX, groupY] = deal(nan(1,max(s.groupId)));
    for i = 1:max(s.groupId)
        groupX(i) = s.groupFcn(x(s.groupId==i));
        groupY(i) = s.groupFcn(y(s.groupId==i));
    end
end


% top histo
subplot(2,2,1); hold on;
set(gca, 'Position', [s.border s.scatPlotSize+s.border s.scatPlotSize 1-s.scatPlotSize-s.border-.01])
if groupNum>1
    for m = 1:groupNum
        kd = ksdensity(x(s.groupId==m), xGrid);
        plot(xGrid, kd, 'LineWidth', 1, 'Color', [s.colors(m,:) .4]);
    end
end
kd = ksdensity(x, xGrid);
plot(xGrid, kd, 'Color', [.2 .2 .2], 'LineWidth', 3);
set(gca, 'XLim', s.xLims, 'Visible', 'off')

% right histo
subplot(2,2,4); hold on;
set(gca, 'Position', [s.scatPlotSize+s.border s.border 1-s.scatPlotSize-s.border-.01 s.scatPlotSize])
if groupNum>1
    for m = 1:groupNum
        kd = ksdensity(y(s.groupId==m), yGrid);
        plot(kd, yGrid, 'LineWidth', 1, 'Color', [s.colors(m,:) .4]);
    end
end
kd = ksdensity(y, yGrid);
plot(kd, yGrid, 'Color', [.2 .2 .2], 'LineWidth', 3);
set(gca, 'YLim', s.yLims, 'Visible', 'off')

% scatter
mainPlot = subplot(2,2,3); hold on;
set(gca, 'Position', [s.border s.border s.scatPlotSize s.scatPlotSize])
if groupNum>1
    for i  = 1:groupNum
        scatter(x(s.groupId==i), y(s.groupId==i), 50, s.colors(i,:), 'filled', 'MarkerFaceAlpha', s.scatAlpha);
    end
else
    scatter(x, y, 50, s.colors, 'filled', 'MarkerFaceAlpha', s.scatAlpha);
end
if ~isempty(s.xlabel); xlabel(s.xlabel); end
if ~isempty(s.xlabel); ylabel(s.ylabel); end
set(gca, 'XLim', s.xLims, 'YLim', s.yLims, 'Box', 'off', 'TickDir', 'out')
uistack(mainPlot, 'bottom')

% add group mean scatters
if groupNum>1; scatter(groupX, groupY, 50, s.colors, 'filled', 'MarkerEdgeColor', 'black'); end


% add crosshair for mean, std for time and distance
if s.showCrossHairs
    if groupNum>1 % average across groups
        line(nanmean(groupX) * [1 1], nanmean(groupY) + [-1 1]*nanstd(groupY), ...
            'linewidth', 5, 'color', 'black')
        line(nanmean(groupX) + [1 -1]*nanstd(groupX), nanmean(groupY)*[1 1], ...
            'linewidth', 5, 'color', 'black')
    else % average across individual samples
        line(nanmean(x) * [1 1], nanmean(y) + [-1 1]*nanstd(y), ...
            'linewidth', 5, 'color', 'black')
        line(nanmean(x) + [1 -1]*nanstd(x), nanmean(y)*[1 1], ...
            'linewidth', 5, 'color', 'black')
    end
end





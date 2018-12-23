function plotAcrossSessions(data, dvs, figTitle)

% settings
addLegend = true;
columns = 1;
mouseSymbols = {'o', '+', '*', '.', 'x', 's', 'd'};

% initializations
rows = ceil(length(dvs) / columns);
figure('name', figTitle, 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 600*columns 200*rows], 'InvertHardcopy', 'off')
mice = unique({data.mouse});
conditions = unique({data.condition});
colors = hsv(length(conditions));
scatters = nan(1,length(mice));


for i = 1:length(mice)
    mouseBins = strcmp({data.mouse}, mice{i});
    for j = 1:length(dvs)
        subplot(rows, columns, j)
        line([data(mouseBins).sessionNum], [data(mouseBins).(dvs{j})], 'color', [.5 .5 .5]); hold on
        [~, colorInds] = ismember({data(mouseBins).condition}, conditions);
        scatter([data(mouseBins).sessionNum], [data(mouseBins).(dvs{j})], ...
            50, colors(colorInds,:), mouseSymbols{i});
    end
end


% pimp figs
for i = 1:length(dvs)
    subplot(rows, columns, i);
    xLims = get(gca, 'xlim');
    set(gca, 'xlim', [0.5 xLims(2)+.5], 'xtick', 1:xLims(2));
    ylabel(dvs{i}, 'interpreter', 'none');
end
xlabel('session number')


if addLegend
    for i = 1:length(mice); scatters(i) = scatter(nan,nan,50,mouseSymbols{i}, 'black'); end % create dummy scatters
    for i = 1:length(conditions); lines(i) = plot([nan nan], 'color', colors(i,:), 'LineWidth', 2); end % create dummy lines
    legend([scatters, lines], cat(2,mice,conditions), 'Location', 'northeastoutside')
end



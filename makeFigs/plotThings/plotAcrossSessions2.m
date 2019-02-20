function plotAcrossSessions2(data, dvs)



% settings
addLegend = true;
columns = 2;
mouseSymbols = {'o', '+', '*', '.', 'x', 's', 'd'};

% initializations
rows = ceil(length(dvs) / columns);
figure('Color', 'white', 'MenuBar', 'none', 'Position', [2000 100 600*columns 150*rows], 'InvertHardcopy', 'off')
mice = unique({data.mouse});
conditions = unique({data.condition});
colors = hsv(length(conditions));
scatters = nan(1,length(mice));
sessions = unique({data.session});
fields = fieldnames(data);


% average across trials for each session
dataTemp = struct();
dataInd = 1;
for i = 1:length(sessions)
    bins = strcmp({data.session}, sessions{i});
    for j = 1:length(fields)
        if ismember(fields{j}, dvs) % copy dv
            dataTemp(dataInd).(fields{j}) = nanmean([data(bins).(fields{j})]);
        else % copy metadata from first row of session
            dataTemp(dataInd).(fields{j}) = [data(find(bins,1,'first')).(fields{j})];
        end
    end
    dataInd = dataInd + 1;
end
data = dataTemp;



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



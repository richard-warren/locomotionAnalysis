function barPlots(data, dvs, figTitle, conditions)

% TO DO: add pairwise significance tests


% settings
addLegend = true;
addStats = true; % peform paired ttests for each pair of conditions
pThresh = .05; % p value threshold for significance
columns = 5;
minConditionNum = 0; % only use a condition after the minConditionNum day of that condition


% initializations
if ~exist('conditions', 'var'); conditions = unique({data.condition}); end
rows = ceil((length(dvs)+addLegend)/columns); % last subplot will contain legend if addLegend
figure('name', figTitle, 'Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 300*columns 250*rows], 'inverthardcopy', 'off')
mice = unique({data.mouse});
xJitters = linspace(-.1,.1,length(mice)); xJitters = xJitters-mean(xJitters); % jitters x position of scatter points
colors = hsv(length(mice));
if addStats
    conditionPairs = nchoosek(1:length(conditions), 2);
    [~, sortInds] = sort(diff(conditionPairs,[],2));
    conditionPairs = conditionPairs(sortInds,:); % sort s.t. more distant comparisons are last
    yPosits = linspace(1.1, 1.1+.025*size(conditionPairs,1), size(conditionPairs,1)); % vertical position of each line, expressed as fraction of y range
end % matrix containing all condition pairs where is row is a pair of conditions


for i = 1:length(dvs)
    
    subplot(rows, columns, i)
    mouseDvAvgs = nan(length(mice), length(conditions)); % last dimension is side of body, only used for dvs broken down into ipsi/contra
    
    for j = 1:length(mice)
        for k = 1:length(conditions)
            bins = strcmp({data.mouse}, mice{j}) & ...
                          strcmp({data.condition}, conditions{k}) & ...
                          [data.conditionNum]>=minConditionNum;
            mouseDvAvgs(j,k) = mean([data(bins).(dvs{i})]);
        end
        
        % plot mouse means    
        line([1:length(conditions)] + xJitters(j), mouseDvAvgs(j,:), 'color', [.2 .2 .2]); hold on
        scatter([1:length(conditions)] + xJitters(j), mouseDvAvgs(j,:), 50, colors(j,:), 'filled', 'MarkerFaceAlpha', .6);
    end
    
    % plot condition means
    for j = 1:length(conditions)
        avg = nanmean(mouseDvAvgs(:,j),1);
        line([j-.1 j+.1], [avg avg], 'linewidth', 3, 'color', 'black')
    end
    
    % add statistics
    if addStats
        yLims = get(gca, 'ylim');
        for j = 1:length(conditionPairs)
            [~,p] = ttest(mouseDvAvgs(:,conditionPairs(j,1)), mouseDvAvgs(:,conditionPairs(j,2)));
            if p<pThresh; lineColor = 'red'; else; lineColor = [.5 .5 .5]; end
            line([conditionPairs(j,:)], yLims(2)*[yPosits(j) yPosits(j)], ...
                'color', lineColor, 'linewidth', 1.0);
        end
        set(gca, 'ylim', [yLims(1) yLims(2)*max(yPosits)]) % set y axis s.t. highest line is at top of y axis
    end
end

% pimp figs
for i = 1:length(dvs)
    subplot(rows, columns, i);
    set(gca, 'xlim', [0.75 length(conditions)+0.25], ...
        'xtick', 1:length(conditions), 'XTickLabel', conditions, 'XTickLabelRotation', 45);
    ylabel(dvs{i}, 'interpreter', 'none');
end

if addLegend
    subplot(rows, columns, length(dvs)+1, 'visible', 'off'); hold on
    for i = 1:length(mice); scatters(i) = scatter(nan,nan,50,colors(i,:), 'filled'); end % create dummy scatters
    legend(scatters, mice);
end

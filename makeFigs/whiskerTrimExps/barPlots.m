function barPlots(data, dvs)

% TO DO: work with ipsi/contra dvs // 


% settings
columns = 3;
minConditionNum = 2; % only use a condition after the minConditionNum day of that condition


% initializations
conditions = unique({data.condition});
rows = ceil(length(dvs)/columns);
figure('Color', 'white', 'MenuBar', 'none', 'Position', [2000 50 300*columns 250*rows], 'inverthardcopy', 'off')
mice = unique({data.mouse});
xJitters = linspace(-.1,.1,length(mice)); xJitters = xJitters-mean(xJitters); % jitters x position of scatter points
colors = hsv(length(mice));
scatters = nan(1,length(mice));


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
        scatters(j) = scatter([1:length(conditions)] + xJitters(j), mouseDvAvgs(j,:), 50, colors(j,:), 'filled');
    end
    
    % plot condition means
    for j = 1:length(conditions)
        avg = nanmean(mouseDvAvgs(:,j),1);
        line([j-.1 j+.1], [avg avg], 'linewidth', 3, 'color', 'black')
        [~,p] = ttest(mouseDvAvgs(:,1), mouseDvAvgs(:,2));
        xlabel(sprintf('p=%.3f', p))
    end
end

% pimp figs
for i = 1:length(dvs)
    subplot(rows, columns, i);
    set(gca, 'xlim', [0.75 length(conditions)+0.25], 'xtick', 1:length(conditions), 'XTickLabel', conditions);
    ylabel(dvs{i});
end

legend(scatters, mice, 'Location', 'northeast', 'color', 'none')

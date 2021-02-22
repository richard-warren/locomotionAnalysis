% explore lick modulations

%% inits

data = getUnitInfo(true);

% add lick importance for every unit (based on additional deviance explained in residual glm models)
data.lickImportance = nan(height(data), 1);
for i = 1:height(data)
    disp(i/height(data))
    fname = fullfile('E:\lab_files\paper2\modelling\glms\residual_glms', ...
        sprintf('%s_cell_%i_glm.mat', data.session{i}, data.unit(i)));
    if exist(fname, 'file')
        load(fname, 'models');
        data.lickImportance(i) = max(0, models{'whiskers', 'dev'});
    end
end

%% for each lick in each SESSION determine: isfirstlick in bout // is reward lick // lick num in bout // time in trial
boutThresh = 1;  % (s) how much time without licking before next lick is considered new bout
rewardPeriod = [0 6];  % (s) time after reward delivery to consider within a reward related lick

sessions = unique(data.session);
inits = repmat({cell(length(sessions),1)}, 1, 4);

lickData = table(inits{:}, 'RowNames', sessions, ...
    'VariableNames', {'isFirstLick', 'isRewardLick', 'boutLickNum', 'lickNum'});

for i = 1:height(lickData)
    fname = fullfile('E:\lab_files\paper2\modelling\runAnalyzed', ...
        sprintf('%s_runAnalyzed.mat', lickData.Properties.RowNames{i}));
    sesData = load(fname, 'rewardTimes', 'lickTimes');
    
    % lick counter
    lickNum = (1:length(sesData.lickTimes))';
    lickData{i, 'lickNum'} = {lickNum};
    
    % whether each lick is first in bout
    timeSinceLastLick = [inf; diff(sesData.lickTimes)];
    isFirstLick = timeSinceLastLick > boutThresh;
    lickData{i, 'isFirstLick'} = {isFirstLick};
    
    % lick number within bout
    boutLickNum = lickNum;
    for j = 1:length(sesData.lickTimes)
        if isFirstLick(j)
            boutLickNum(j:end) = boutLickNum(j:end) - boutLickNum(j) + 1;
        end
    end
    lickData{i, 'boutLickNum'} = {boutLickNum};
    
    % whether each lick is within a reward period
    isRewardLick = any((sesData.lickTimes > (sesData.rewardTimes+rewardPeriod(1))') & ...
                       (sesData.lickTimes < (sesData.rewardTimes+rewardPeriod(2))'), 2);
    lickData{i, 'isRewardLick'} = {isRewardLick};
end


%% psth grids (lick resp X lick num for reward and non-reward licks // lick resp X time in trial)


% settings
ncols = 8;  % divisible by length(paws) ideally
tbins = 3;  % divide trials into this many bins
plots = {'vel_response', 'lick_response', 'response', 'predicted_full', 'predicted_noreward', 'lick_psth'};
xlims = [x(1) x(end)];
figpos = [2.00 475.00 2554.00 881.00];
colors = lines(length(plots));
ylimMatchRows = [3,4];  % match y limits for these rows within each column

% inits
mask = linspace(1, 0, tbins);
close all
figure('name', 'reward responses', 'color', 'white', 'position', figpos, 'menubar', 'none'); hold on
fignum = 1; offset = 0;


tic
for i = 1:height(data)
    
    if (i-offset) > ncols
        offset = offset + ncols;
        saveas(gcf, ['E:\lab_files\paper2\plots\rewards\overtime_psths\' num2str(fignum) '.png'])
        
        set(gcf, 'PaperOrientation', 'landscape');
        print(gcf, ['E:\lab_files\paper2\plots\rewards\psths' num2str(fignum) '.pdf'], '-dpdf', '-bestfit')
        
        fignum = fignum + 1;
        figure('name', 'reward responses', 'color', 'white', 'position', figpos); hold on
    end
    
    subplotNum = i-offset;
    
    ylims = nan(length(plots),2);
    
    for j = 1:length(plots)
        subplot(length(plots), ncols, subplotNum + (j-1)*ncols); hold on
        resp = data{i, plots{j}}{1};
        
        % get bins
        slice = resp(:,1);  % slice of firing rate first time bin across trails
        startInd = find(~isnan(slice), 1, 'first');
        endInd = find(~isnan(slice), 1, 'last');
        binNums = nan(1, length(slice));
        binNums(startInd:endInd) = round(linspace(1, tbins, sum(~isnan(slice))));
        
        for k = 1:tbins
            if ~isempty(resp)
                bins = k==binNums;
                mn = nanmean(resp(bins,:), 1);
                stdev = nanstd(resp(bins,:), 1);
                patch([x fliplr(x)], [(-stdev+mn) fliplr(stdev+mn)], colors(j,:).*mask(k), ...
                    'FaceAlpha', .25, 'EdgeColor', 'none')           % shaded error bars
                plot(x, mn, 'color', colors(j,:).*mask(k), 'lineWidth', 2)  % mean
                plot([0 0], ylim, 'color', [0 0 0 .4])               % vertical line at t=0
            end
        end
        
%         yticks = get(gca, 'ytick');
        set(gca, 'XLim', xlims, 'xcolor', 'none')  % 'YTick', yticks([1,end]), 
        ylims(j,:) = ylim;
        
        if j==1
            tit1 = sprintf('%s(%i)', data{i, 'session'}{1}, data{i, 'unit'});
            tit2 = data{i, 'nucleus'}{1};
            title({tit1, tit2}, ...
                'Interpreter', 'none', 'FontWeight', 'normal', 'FontSize', 8)
        end
        
        if subplotNum==1
            ylabel(plots{j}, 'Interpreter', 'none', 'FontWeight', 'bold')
            if j==length(plots); set(gca, 'xcolor', get(gca, 'ycolor')); end
        end
    end
    
    ylimMatched = [min(ylims(ylimMatchRows,1)), max(ylims(ylimMatchRows,2))];
    for j = ylimMatchRows
        subplot(length(plots), ncols, subplotNum + (j-1)*ncols); hold on
        set(gca, 'ylim', ylimMatched)
    end
    
end

set(gca, 'xcolor', get(gca, 'ycolor'))  % time axis visible only for final subplot
saveas(gcf, ['E:\lab_files\paper2\plots\rewards\overtime_psths\' num2str(fignum) '.png'])
toc
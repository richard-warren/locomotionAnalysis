% explore lick modulations

% idea: should look at first lick after reward when anticipatory licks were
% already going // this could test whether signal really necessary to 'kick
% off' the cpg...

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
        data.lickImportance(i) = max(0, models{'lick', 'dev'});
    end
end

% for each lick in each SESSION determine: isfirstlick in bout // is reward lick // lick num in bout // time in trial
boutThresh = 1;  % (s) how much time without licking before next lick is considered new bout
rewardPeriod = [0 6];  % (s) time after reward delivery to consider within a reward related lick

sessions = unique(data.session);
inits = repmat({cell(length(sessions),1)}, 1, 6);

lickData = table(inits{:}, 'RowNames', sessions, ...
    'VariableNames', {'time', 'isFirstLick', 'isRewardLick', 'boutLickNum', 'lickNum', 'hasAntipLicks'});

for i = 1:height(lickData)
    fname = fullfile('E:\lab_files\paper2\modelling\runAnalyzed', ...
        sprintf('%s_runAnalyzed.mat', lickData.Properties.RowNames{i}));
    sesData = load(fname, 'rewardTimes', 'lickTimes');
    lickData{i, 'time'} = {sesData.lickTimes};
    
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
    
    % for each reward lick, code whether anticipatory licks were present for that trial
    lickData{i, 'hasAntipLicks'} = {false(length(sesData.lickTimes),1)};
    for j = 1:length(sesData.rewardTimes)
        if j<length(sesData.rewardTimes); maxt = sesData.rewardTimes(j+1); else; maxt = inf; end
        lickBins = (sesData.lickTimes > sesData.rewardTimes(j)) & ...
                   (sesData.lickTimes < maxt) & ...
                   isRewardLick;
        hasAntipLicks = any((sesData.lickTimes > (sesData.rewardTimes(j)-boutThresh)) & ...
                            (sesData.lickTimes < (sesData.rewardTimes(j))));
        lickData{i, 'hasAntipLicks'}{1}(lickBins) = hasAntipLicks;
    end
end


%% psth grids (lick resp X lick num for reward and non-reward licks // lick resp X time in trial)

% todo:
% make sure only spans recording period for each unit
% check against vid i'm sure of
% subtract mean to handle offsets

% settings
ncols = 8;  % divisible by length(paws) ideally
x = linspace(-.1, .1, 200);
figpos = [2.00 717.00 2551.00 693.00];
maxLickNum = 4;
numBins = 3;  % bins for time in trial
colors = colorme(3, 'saturation', .8, 'offset', .25, 'showSamples', false);
dcRemove = true;  % whether to subtract mean of every response to remove dc shifts

% inits
nrows = 3;
timeMask = linspace(1, 0, numBins);
numMask = linspace(1, 0, maxLickNum);
close all
figure('name', 'reward responses', 'color', 'white', 'position', figpos, 'menubar', 'none'); hold on
fignum = 1; offset = 0;
prevSession = '';
xlims = [x(1) x(end)];

% curate data
datasub = data(ismember(data.nucleus, {'fastigial', 'interpositus', 'dentate'}),:);
[~, sortInds] = sortrows(datasub(:, {'lickImportance', 'nucleus'}));
sortInds = flipud(sortInds);
datasub = datasub(sortInds,:);

tic
for i = 1:height(datasub)
    
    % make new fig if necessary
    if (i-offset) > ncols
        offset = offset + ncols;
        saveas(gcf, ['E:\lab_files\paper2\plots\licks\psths' num2str(fignum) '.png'])
        fignum = fignum + 1;
        figure('name', 'reward responses', 'color', 'white', 'position', figpos, 'menubar', 'none'); hold on
    end
    subplotNum = i-offset;
    
    % load new session data if necessary
    if ~strcmp(prevSession, datasub.session{i})
        neuralData = load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'neuralData', [datasub.session{i} '_neuralData.mat']), ...
            'unit_ids', 'spkRates', 'timeStamps');
        sesInd = find(strcmp(lickData.Properties.RowNames, datasub.session{i}));  % index for lickData, which has one row per session
        prevSession = datasub.session{i};
    end
    
    % get psth for unit
    unitInd = find(datasub.unit(i) == neuralData.unit_ids);
    resp = interp1(neuralData.timeStamps, neuralData.spkRates(unitInd,:), lickData.time{sesInd} + x);  % one row per trial
    ylims = nan(3,2);  % store ylims for all plots so we can match them at the end
    
    % find bins for licks occuring during recording for unit
    tsub = neuralData.timeStamps(~isnan(neuralData.spkRates(unitInd,:)));  % times when unit is active
    unitLickBins = lickData.time{sesInd}>tsub(1) & lickData.time{sesInd}<tsub(end);
    
    
    % 1) bout lick num for rewarded licks (only reward licks without anticipatory licks
    subplot(nrows, ncols, subplotNum); hold on
    bins = lickData{sesInd, 'isRewardLick'}{1} & ~lickData{sesInd, 'hasAntipLicks'}{1};  % licks that are in reward period for trials in which there are not anticipatory licks
    groups = lickData{sesInd, 'boutLickNum'}{1} .* bins;  % set invalid licks to zero
        
    plts = nan(1, maxLickNum);
    for j = 1:maxLickNum
        bins = groups==j;
        mn = nanmean(resp(bins,:), 1);
        if dcRemove; mn = mn - mean(mn); end
        st = nanstd(resp(bins,:), 1) / sqrt(sum(bins));
        patch([x fliplr(x)], [(-st+mn) fliplr(st+mn)], colors(1,:).*numMask(j), ...
            'FaceAlpha', .25, 'EdgeColor', 'none')           % shaded error bars
        plts(j) = plot(x, mn, 'color', colors(1,:).*numMask(j), 'lineWidth', 2);  % mean
    end
    
    ylims(1,:) = ylim;
    if subplotNum==1
        legend(plts, split(num2str(1:maxLickNum)), 'autoupdate', 'off', 'fontsize', 8)
        ylabel('bout num (reward licks)');
    end
    
    tit1 = sprintf('%s(%i)', datasub{i, 'session'}{1}, datasub{i, 'unit'});
    tit2 = datasub{i, 'nucleus'}{1};
    title({tit1, tit2}, ...
        'Interpreter', 'none', 'FontWeight', 'normal', 'FontSize', 8)
    
    
    
    % 2) bout lick num for non reward lick (only non reward licks)
    subplot(nrows, ncols, subplotNum + ncols); hold on
    bins = ~lickData{sesInd, 'isRewardLick'}{1};  % licks that are in reward period for trials in which there are not anticipatory licks
    groups = lickData{sesInd, 'boutLickNum'}{1} .* bins;  % set invalid licks to zero
        
    plts = nan(1, maxLickNum);
    for j = 1:maxLickNum
        bins = groups==j;
        mn = nanmean(resp(bins,:), 1);
        if dcRemove; mn = mn - mean(mn); end
        st = nanstd(resp(bins,:), 1) / sqrt(sum(bins));
        patch([x fliplr(x)], [(-st+mn) fliplr(st+mn)], colors(2,:).*numMask(j), ...
            'FaceAlpha', .25, 'EdgeColor', 'none')           % shaded error bars
        plts(j) = plot(x, mn, 'color', colors(2,:).*numMask(j), 'lineWidth', 2);  % mean
    end
    ylims(2,:) = ylim;
    if subplotNum==1
        legend(plts, split(num2str(1:maxLickNum)), 'autoupdate', 'off', 'fontsize', 8)
        ylabel('bout num (antic licks)');
    end
    
    
    % 3) lick resp by time in session (only reward licks;
    subplot(nrows, ncols, subplotNum + 2*ncols); hold on
    bins = lickData{sesInd, 'isRewardLick'}{1} & unitLickBins;
    lickNum = lickData.lickNum{sesInd};
    lickNum(~bins) = nan;
    [~, ~, groups] = histcounts(lickNum, numBins);
    
    plts = nan(1, numBins);
    for j = 1:numBins
        bins = groups==j;
        mn = nanmean(resp(bins,:), 1);
        if dcRemove; mn = mn - mean(mn); end
        st = nanstd(resp(bins,:), 1) / sqrt(sum(bins));
        patch([x fliplr(x)], [(-st+mn) fliplr(st+mn)], colors(3,:).*timeMask(j), ...
            'FaceAlpha', .25, 'EdgeColor', 'none')           % shaded error bars
        plts(j) = plot(x, mn, 'color', colors(3,:).*timeMask(j), 'lineWidth', 2);  % mean
    end
    ylims(3,:) = ylim;
    if subplotNum==1
        legend(plts, [{'early'}, repmat({''},1,numBins-2), {'late'}], 'autoupdate', 'off', 'fontsize', 8)
        ylabel('time in session');
    end
    
    % set axis props
    ylims = [min(ylims(:,1)), max(ylims(:,2))];
    for j = 1:3
        subplot(nrows, ncols, subplotNum + (j-1)*ncols); hold on
        plot([0 0], ylims, 'color', [0 0 0 .4])
        set(gca, 'xlim', xlims, 'ylim', ylims)
        limitticks
    end
end

set(gca, 'xcolor', get(gca, 'ycolor'))  % time axis visible only for final subplot
saveas(gcf, ['E:\lab_files\paper2\plots\licks\psths' num2str(fignum) '.png'])
toc
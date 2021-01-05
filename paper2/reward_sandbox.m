%% explore reward responses

%% collect data
% tbl with real responses to reward, and responses predicted by full models
% and models without reward predictor


data = getUnitInfo();
sessions = unique(data.session);
n = height(data);


% settings
x = linspace(-.5, 2, 200);  % (s) x grid for responses
init = repmat({nan(1,length(x))}, n, 1);
tbl = table(init, init, init, init, init, nan(n,1), nan(n,1), nan(n,1), nan(n,1), 'VariableNames', ...
    {'response', 'predicted_full', 'predicted_noreward', 'lick_response', 'vel_response', 'dev_noreward', 'dev_full', 'mean', 'std'});

for i = 1:length(sessions)
    disp(i/length(sessions))
    % load ephys info
    neuralData = load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'neuralData', [sessions{i} '_neuralData.mat']), ...
        'unit_ids', 'spkRates', 'timeStamps');
    unit_ids = neuralData.unit_ids;
    load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'predictors', [sessions{i} '_predictors.mat']), 'predictors');
    rewardTimes = predictors{'reward_normal', 'data'}{1};
    X = rewardTimes + x;  % each row is the time grid for a single trial
    
    % lick rate
    licks = predictors{'lick', 'data'}{1};
    [lickRate, lickTimes] = getFiringRate(licks, 'kernel', 'gauss', 'kernelSig', .05);
    lickResponses = interp1(lickTimes, lickRate, X);
    
    % velocity
    velResponses = interp1(predictors{'velocity', 't'}{1}, predictors{'velocity', 'data'}{1}, X);
        
    for j = 1:length(unit_ids)
        
        % load models
        rowInd = find(strcmp(data.session, sessions{i}) & data.unit==unit_ids(j));
        fname = ['E:\lab_files\paper2\modelling\glms\upper_lower_glms\' sessions{i} '_cell_' num2str(unit_ids(j)) '_glm.mat'];
        if exist(fname, 'file')
            models = load(fname);
            t = models.fitdata.t;
            dt = t(2) - t(1);

            % neural response
            tbl{rowInd, 'response'}{1} = interp1(neuralData.timeStamps, neuralData.spkRates(j,:), X);
            tbl{rowInd, 'mean'} = nanmean(neuralData.spkRates(j,:));
            tbl{rowInd, 'std'} = nanmean(neuralData.spkRates(j,:));

            % predicted response (full)
            yhat = exp(models.models{'full', 'model_in'}{1}.fit_preval) / dt;
            tbl{rowInd, 'predicted_full'}{1} = interp1(t, yhat, X);

            % predicted response (no reward)
            yhat = exp(models.models{'reward', 'model_out'}{1}.fit_preval) / dt;
            tbl{rowInd, 'predicted_noreward'}{1} = interp1(t, yhat, X);

            % deviance explained
            tbl{rowInd, 'dev_noreward'} = models.models{'full', 'dev_in'} - models.models{'reward', 'dev_out'};
            tbl{rowInd, 'dev_full'} = models.models{'full', 'dev_in'};
            
            % lick rates and vel
            tbl{rowInd, 'lick_response'}{1} = lickResponses;  % todo: should be done with one linear across all units for session, or saved once pre session to avoid redundancy
            tbl{rowInd, 'vel_response'}{1} = velResponses;
        end
    end
end
data = cat(2, data, tbl);
disp('all done!')

%% psth grids

% settings
ncols = 8;  % divisible by length(paws) ideally
nrows = 6;
colors = [.4 .4 .4; lines(2)];  % response // predicted full // predicted no reward
xlims = [x(1) x(end)];


% inits
% nrows = ceil(height(data) / ncols);
close all
figure('name', 'reward responses', 'color', 'white', 'position', [1.00 1.00 2560.00 1363.00]); hold on
plots = {'response', 'predicted_full', 'predicted_noreward'};
fignum = 1; offset = 0;


tic
for i = 1:10%:height(data)
    
    
    if (i-offset) > ncols*nrows
        offset = offset + ncols*nrows;
        saveas(gcf, ['E:\lab_files\paper2\plots\rewards\psths' num2str(fignum) '.png'])
        
        set(gcf, 'PaperOrientation', 'landscape');
        print(gcf, ['E:\lab_files\paper2\plots\rewards\psths' num2str(fignum) '.pdf'], '-dpdf', '-bestfit')
        
        fignum = fignum + 1;
        figure('name', 'reward responses', 'color', 'white', 'position', [1.00 1.00 2560.00 1363.00]); hold on
    end
    
    subplot(nrows, ncols, i-offset); hold on
    
    for j = 1:length(plots)
        resp = data{i, plots{j}}{1};
        if ~isempty(resp)
            mn = nanmean(resp, 1);
            stdev = nanstd(resp, 1);
            patch([x fliplr(x)], [(-stdev+mn) fliplr(stdev+mn)], colors(j,:), ...
                'FaceAlpha', .25, 'EdgeColor', 'none')           % shaded error bars
            plot(x, mn, 'color', colors(j,:), 'lineWidth', 2)  % mean
            plot([0 0], ylim, 'color', [0 0 0 .4])               % vertical line at t=0
        end
    end
    
    tit1 = sprintf('%s(%i)', data{i, 'session'}{1}, data{i, 'unit'});
    tit2 = data{i, 'nucleus'}{1};
    title({tit1, tit2}, ...
        'Interpreter', 'none', 'FontWeight', 'normal', 'FontSize', 8)
    yticks = get(gca, 'ytick');
    set(gca, 'XLim', xlims, 'YTick', yticks([1,end]), 'xcolor', 'none')
end

set(gca, 'xcolor', get(gca, 'ycolor'))  % time axis visible only for final subplot
saveas(gcf, ['E:\lab_files\paper2\plots\rewards\psths' num2str(fignum) '.png'])

set(gcf, 'PaperOrientation', 'landscape'); print(gcf, ['E:\lab_files\paper2\plots\rewards\psths' num2str(fignum) '.pdf'], '-dpdf', '-bestfit')
toc


%% scatters of dev with and without rewards, colored by nucleus

close all;
figure('color', 'white', 'position', [124.00 882.00 1138.00 316.00], 'menubar', 'none');

subplot(1,3,1); hold on
scatter(data{:, 'dev_noreward'}, data{:, 'dev_full'}, [], 'black', 'filled', 'MarkerFaceAlpha', .4)
plot(xlim, xlim, 'Color', [0 0 0 .4])
xlabel('no reward'); ylabel('full')

subplot(1,3,2)
mat = data{:, {'dev_noreward', 'dev_full'}}';
deltas = diff(mat, [], 1);
histogram(deltas, 'Normalization', 'probability', 'FaceColor', 'black', 'EdgeColor', 'none')
xlabel('\Delta deviance explained')
set(gca, 'box', 'off')

subplot(1,3,3)
barFancy(mat, 'ylabel', 'deviance explained', 'levelNames', {{'no reward', 'full'}}, ...
    'comparisons', [1 2], 'showViolins', true, 'showBars', false, 'barWidth', .5, ...
    'scatterColors', [0 0 0], 'scatterAlpha', .05, 'showErrorBars', false)

saveas(gcf, 'E:\lab_files\paper2\plots\rewards\scatters.png')


%% heatmaps of rewards, noreward predictions, and residuals!

plots = {'response', 'predicted_noreward', 'residual'};
close all
% cmap = customcolormap([0 .5 1], [1 .2 .2; 1 1 1; .2 .2 1]);
figure('color', 'white', 'position', [2.00 722.00 1278.00 688.00], 'menubar', 'none');

resp = cellfun(@(x) nanmean(x,1), data{:, 'response'}, 'UniformOutput', false); resp = cat(1, resp{:});
pred_noreward = cellfun(@(x) nanmean(x,1), data{:, 'predicted_noreward'}, 'UniformOutput', false); pred_noreward = cat(1, pred_noreward{:});
pred_residual = resp - pred_noreward;
responses = {resp, pred_noreward, pred_residual};
bins = ~isnan(resp(:,1));
zeroind = knnsearch(x', 0);

% get sort inds
[~, maxInds] = max(resp, [], 2);
[~, sortInds] = sort(maxInds);

for i = 1:length(responses)
    subplot(1, length(plots), i); hold on
    title(plots{i}, 'Interpreter', 'none')
    hmap = responses{i};
%     hmap = (hmap - hmap(:,zeroind));
%     hmap = (hmap - responses{1}(:,1)) ./ data{:, 'std'};
    
    imagesc(x, 1:sum(bins), hmap(sortInds(bins),:))
    if i==1; clim = get(gca, 'clim'); end
%     colormap(gca, 'bone')
    plot([0 0], [1 sum(bins)], 'color', 'black')
    set(gca, 'box', 'off', 'ycolor', 'none', 'xlim', [x(1) x(end)], 'CLim', clim)
end

saveas(gcf, 'E:\lab_files\paper2\plots\rewards\heatmaps.png')


%% psth grids

% settings
ncols = 8;  % divisible by length(paws) ideally
tbins = 3;  % divide trials into this many bins
plots = {'vel_response', 'lick_response', 'response', 'predicted_noreward'};
xlims = [x(1) x(end)];
figpos = [2.00 821.00 2554.00 535.00];
colors = lines(length(plots));


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
    
    slice = data{i, 'response'}{1}(:,1);  % slice of firing rate first time bin across trails
    startInd = find(~isnan(slice), 1, 'first');
    endInd = find(~isnan(slice), 1, 'last');
    binNums = nan(1, length(slice));
    binNums(startInd:endInd) = round(linspace(1, tbins, sum(~isnan(slice))));
    
    for j = 1:length(plots)
        subplot(length(plots), ncols, subplotNum + (j-1)*ncols); hold on
        resp = data{i, plots{j}}{1};
        
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
        
        yticks = get(gca, 'ytick');
        set(gca, 'XLim', xlims, 'YTick', yticks([1,end]), 'xcolor', 'none')
        
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
    
end

set(gca, 'xcolor', get(gca, 'ycolor'))  % time axis visible only for final subplot
saveas(gcf, ['E:\lab_files\paper2\plots\rewards\overtime_psths\' num2str(fignum) '.png'])
toc

%% find matched slowdowns!

% settings
session = sessions{50};
xlims = [-2 1];
nbest = 20;  % find nbest best matches
c = lines(2); c = c(2,:);

% inits
load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'predictors', [session '_predictors.mat']), 'predictors');
vel = predictors{'velocity', 'data'}{1};
t = predictors{'velocity', 't'}{1};
rewardTimes = predictors{'reward_all', 'data'}{1};

% get average response
dt = nanmedian(diff(t));
x = xlims(1) : dt : xlims(2);
X = x + rewardTimes;
responses = interp1(t, vel, X);
kernel = nanmean(responses, 1);

% mask out peri-event periods
maskInds = knnsearch(t', X(:));
vel(maskInds) = nan;


n = length(kernel);
diffs = nan(1, length(vel));
for i = 1 : length(vel)-n
%     diffs(i) = mean(abs(vel(i:i+n-1) - kernel));
    diffs(i) = sum((vel(i:i+n-1) - kernel).^2);
end

figure('color', 'white', 'position', [41.00 812.00 1186.00 411.00]);

% mean slow down
subplot(2,3,1); hold on
plot(x, responses', 'color', [0 0 0 .1])
plot(x, mean(responses,1), 'LineWidth', 3, 'color', [0 0 0])
set(gca, 'xlim', xlims);
title('kernel')

% raw signal
subplot(2,3,4:6); hold on
ylims = ylim;
plot(repmat(rewardTimes,1,2), ylim, 'color', [.2 .2 .8])  % plot reward times
plot(t, vel, 'color', [.4 .4 .4])
set(gca, 'ylim', ylims)

% find and plot peaks
[~, peaks] = findpeaks(-diffs, t, 'SortStr', 'descend', 'MinPeakDistance', 5, 'NPeaks', nbest);
peaks = peaks - x(1);  % shift to the left to compensate for how diffs was computed from left edge of kernel
for i = 1:length(peaks)
    xsub = x + peaks(i);
    plot(xsub, kernel, 'color', c, 'linewidth', 2)
end

% matched mean slow down
subplot(2,3,2); hold on
matched = interp1(t, vel, x + peaks');
plot(x, matched', 'color', [c .1])
plot(x, mean(matched, 1), 'LineWidth', 3, 'color', c)
set(gca, 'xlim', xlims);
title('matched')

% matched over kernel
subplot(2,3,3); hold on

mn = mean(responses, 1);
plot(x, mn, 'LineWidth', 3, 'color', [0 0 0])
stdev = std(responses, 1);
patch([x fliplr(x)], [(-stdev+mn) fliplr(stdev+mn)], [0 0 0], 'FaceAlpha', .25, 'EdgeColor', 'none')

mn = mean(matched, 1);
plot(x, mn, 'LineWidth', 3, 'color', c)
stdev = std(matched, 1);
patch([x fliplr(x)], [(-stdev+mn) fliplr(stdev+mn)], c, 'FaceAlpha', .25, 'EdgeColor', 'none')


set(gca, 'xlim', xlims);
title('kernel + matched')










































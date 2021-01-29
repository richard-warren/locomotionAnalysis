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
tbl = table(init, init, init, init, init, init, ...
    nan(n,1), nan(n,1), nan(n,1), nan(n,1), 'VariableNames', ...
    {'response', 'predicted_full', 'predicted_noreward', 'lick_response', 'vel_response', 'lick_psth', ...
    'dev_noreward', 'dev_full', 'mean', 'std'});

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
    Xlick = licks + x;
    
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
            
            % lick psth
            tbl{rowInd, 'lick_psth'}{1} = interp1(neuralData.timeStamps, neuralData.spkRates(j,:), Xlick);
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


%% psth grids (vel, lick reward, reward_predictions, lick psth)

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

%% matched lick bouts and slow downs!

% settings
velWindow = [-2 6];
lickWindow = [0 4];
nbestVel = 20;
nbestLick = 20;
colors = flipud(lines(2));
ncols = 10;
figpos = [2.00 563.00 2558.00 793.00];
corrMaxLag = 2;  % (s)

% inits
% data = getUnitInfo();
xLims = [min(velWindow(1), lickWindow(1)), max(velWindow(2), lickWindow(2))];
x = linspace(xLims(1), xLims(2), 1000);
corrDt = .01;
maxLagSmps = corrMaxLag/corrDt;

close all
figure('name', 'reward responses', 'color', 'white', 'position', figpos); hold on
previousSession = '';
figInd = 1;

% for session
for i = 1:height(data)
    session = data{i, 'session'}{1};
    if ~strcmp(previousSession, session)
        neuralData = load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'neuralData', [session '_neuralData.mat']));
        
        load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'predictors', [session '_predictors.mat']), 'predictors');
        rewardNormal = predictors{'reward_normal', 'data'}{1};
        rewardSurprise = predictors{'reward_surprise', 'data'}{1};
        rewardOmission = predictors{'reward_omission', 'data'}{1};
        
        conditions = {'normal', 'surprise', 'omission'};
        rewardTimes = {rewardNormal, rewardSurprise, rewardOmission};
        rewardColors = [0 0 0 1; 0 .6 0 .5; .8 0 0 .5];
        
        lickTimes = predictors{'lick', 'data'}{1};
        [lickRate, lickRateTimes] = getFiringRate(lickTimes, 'kernel', 'gauss', 'kernelSig', .05, ...
            'fs', 100, 'tLims', [neuralData.timeStamps(1) neuralData.timeStamps(end)]);
        vel = predictors{'velocity', 'data'}{1};
        velTimes = predictors{'velocity', 't'}{1};
        
        velMatches = findMatchedSignalTimes(vel, velTimes, rewardNormal, velWindow, nbestVel);
        lickMatches = findMatchedSignalTimes(lickRate, lickRateTimes, rewardNormal, lickWindow, nbestLick);
    end
    unitInd = find(neuralData.unit_ids==data{i, 'unit'});
    colInd = mod(i-1, ncols)+1;
    
    % compute smoother firing rates for matches epochs, bc far fewer trials
    if any(neuralData.spkTimes{unitInd})
        [smoothedSpks, smoothedSpksTimes] = getFiringRate(...
            neuralData.spkTimes{unitInd}, 'kernel', 'gauss', 'kernelSig', .06);

        % vel
        subplot(5,ncols,colInd)
        resp = interp1(velTimes, vel, rewardNormal + x);
        shadedErrorBar(x, resp, {@nanmean, @nanstd}, 'lineProps', {'lineWidth', 3, 'color', [0 0 0]})
        resp = interp1(velTimes, vel, velMatches + x);
        shadedErrorBar(x, resp, {@nanmean, @nanstd}, 'lineProps', {'lineWidth', 3, 'color', colors(1,:)})
        if colInd==1; ylabel('velocity'); end

        tit1 = sprintf('%s(%i)', data{i, 'session'}{1}, data{i, 'unit'});
        tit2 = data{i, 'nucleus'}{1};
        title({tit1, tit2}, 'Interpreter', 'none', 'FontWeight', 'normal', 'FontSize', 8)

        % lick rate
        subplot(5,ncols,colInd+ncols)
        resp = interp1(lickRateTimes, lickRate, rewardNormal + x);
        shadedErrorBar(x, resp, {@nanmean, @nanstd}, 'lineProps', {'lineWidth', 3, 'color', [0 0 0]})
        resp = interp1(lickRateTimes, lickRate, lickMatches + x);
        shadedErrorBar(x, resp, {@nanmean, @nanstd}, 'lineProps', {'lineWidth', 3, 'color', colors(2,:)})
        if colInd==1; ylabel('lick rate'); end

        % firing rate
        subplot(5,ncols,colInd+2*ncols); hold on
        resp = interp1(neuralData.timeStamps, neuralData.spkRates(unitInd,:), rewardNormal + x);
        shadedErrorBar(x, resp, {@nanmean, @nanstd}, 'lineProps', {'lineWidth', 2, 'color', [0 0 0]})
        resp = interp1(smoothedSpksTimes, smoothedSpks, velMatches + x);
        plot(x, nanmean(resp,1), 'color', [colors(1,:) .75], 'LineWidth', 2)
        resp = interp1(smoothedSpksTimes, smoothedSpks, lickMatches + x);
        plot(x, nanmean(resp,1), 'color', [colors(2,:) .75], 'LineWidth', 2)
        if colInd==1; ylabel('firing rate'); xlabel('time after reward (s)'); end
        
        % xcorr
        subplot(5,ncols,colInd+3*ncols); hold on
        plot(corrMaxLag*[-1 1], [0 0], 'color', [1 1 1]*.6)  % line at y=0
        plot([0 0], [-1 1], 'color', [1 1 1]*.6)            % line at x=0
        
        t = neuralData.timeStamps(~isnan(neuralData.spkRates(unitInd,:)));
        t = t(1) : corrDt : t(end);  % interpole onto common .01 s time grid
        velInterp = zscore(interp1(velTimes, vel, t));
        lickRateInterp = zscore(interp1(lickRateTimes, lickRate, t));
        frInterp = zscore(interp1(neuralData.timeStamps, neuralData.spkRates(unitInd,:), t));
        
        [corrVel, lags] = xcorr(frInterp, velInterp, maxLagSmps, 'unbiased');
        corrLick = xcorr(frInterp, lickRateInterp, maxLagSmps, 'unbiased');
        plot(lags*corrDt, corrVel, 'color', colors(1,:), 'LineWidth', 2)
        plot(lags*corrDt, corrLick, 'color', colors(2,:), 'LineWidth', 2)
        set(gca, 'ylim', [-1 1])
        if colInd==1; ylabel('correlation'); xlabel('\leftarrow brain leads (lag [s]) body leads \rightarrow'); end
        
        % fr responses to different reward types
        subplot(5,ncols,colInd+4*ncols); hold on
        for j = 1:3
            resp = interp1(smoothedSpksTimes, smoothedSpks, rewardTimes{j} + x);
            shadedErrorBar(x, resp, {@nanmean, @nanstd}, 'patchSaturation', .1, ...
                'lineProps', {'lineWidth', 3, 'color', rewardColors(j,:)})
        end
        if colInd==1; ylabel('firing rate'); xlabel('time after reward (s)'); end
        
    else
        fprintf('WARNING! No spikes for %s unit %i\n', session, data{i, 'unit'})
    end
    
    % make new figure if necessary
    if mod(i, ncols)==0 || i==height(data)
        saveas(gcf, ['E:\lab_files\paper2\plots\rewards\reward_matching\' num2str(figInd) '.png'])
        if i<height(data)
            figure('name', 'reward responses', 'color', 'white', 'position', figpos); hold on
            figInd = figInd + 1;
        end
    end
end

%% responses to different reward types (fr, vel, lick rate)


% settings
xLims = [-2 6];
ncols = 10;  % 10
colors = [0 0 0 1; 0 .6 0 .5; .8 0 0 .5];
figpos = [2.00 827.00 2558.00 529.00];
sortCols = {'nucleus', 'session'};


% inits
% data = getUnitInfo();  % only run once
x = linspace(xLims(1), xLims(2), 1000);
[~, sortInds] = sortrows(data, sortCols);

close all
figure('name', 'reward responses', 'color', 'white', 'position', figpos); hold on
previousSession = '';
figInd = 1;

% for session
for i = 1:height(data)
    session = data{sortInds(i), 'session'}{1};
    if ~strcmp(previousSession, session)
        neuralData = load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'neuralData', [session '_neuralData.mat']));
        
        load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'predictors', [session '_predictors.mat']), 'predictors');
        rewardNormal = predictors{'reward_normal', 'data'}{1};
        rewardSurprise = predictors{'reward_surprise', 'data'}{1};
        rewardOmission = predictors{'reward_omission', 'data'}{1};
        
        conditions = {'normal', 'surprise', 'omission'};
        rewardTimes = {rewardNormal, rewardSurprise, rewardOmission};
        
        
        lickTimes = predictors{'lick', 'data'}{1};
        [lickRate, lickRateTimes] = getFiringRate(lickTimes, 'kernel', 'gauss', 'kernelSig', .05, ...
            'fs', 100, 'tLims', [neuralData.timeStamps(1) neuralData.timeStamps(end)]);
        vel = predictors{'velocity', 'data'}{1};
        velTimes = predictors{'velocity', 't'}{1};
    end
    
    unitInd = find(neuralData.unit_ids==data{sortInds(i), 'unit'});
    colInd = mod(i-1, ncols)+1;
    
    if any(neuralData.spkTimes{unitInd})
        
        % compute smoother firing rates for matches epochs, bc far fewer trials
        [smoothedSpks, smoothedSpksTimes] = getFiringRate(...
            neuralData.spkTimes{unitInd}, 'kernel', 'gauss', 'kernelSig', .06);
        
        % times all times and firing rate into tableR
        unitData = table({velTimes, lickRateTimes, smoothedSpksTimes}', ...
                         {vel, lickRate, smoothedSpks}', ...
                         'VariableNames', {'t', 'data'}, ...
                         'RowNames', {'velocity', 'lick rate', 'firing rate'});
        
        for row = 1:height(unitData)
            subplot(height(unitData), ncols, colInd+(row-1)*ncols); hold on
            
            for j = 1:3  % for each of 3 response types
                rewardBins = rewardTimes{j}>=smoothedSpksTimes(1) & rewardTimes{j}<=smoothedSpksTimes(end);  % bins where unit is not nan
                resp = interp1(unitData{row, 't'}{1}, ...
                               unitData{row, 'data'}{1}, ...
                               rewardTimes{j}(rewardBins) + x);
                mn = nanmean(resp,1);
                dev = nanstd(resp,1);
                patch([x fliplr(x)], [mn+dev fliplr(mn-dev)], colors(j,1:3), ...
                    'EdgeColor', 'none', 'FaceAlpha', .1)
                plot(x, mn, 'Color', colors(j,:), 'LineWidth', 3)
            end
            yLims = ylim;
            plot([0 0], yLims, 'Color', [0 0 0 .5])  % line at x=0
            set(gca, 'ylim', yLims)
            
            if colInd==1
                ylabel(unitData.Properties.RowNames{row});
                xlabel('time after reward (s)');
                if row==height(unitData)
                    for k = 1:3; lines(k) = plot([nan nan], 'color', colors(k,:), 'LineWidth', 2); end % create dummy lines
                    legend(lines, conditions, 'Location', 'best', 'Box', 'off', ...
                        'AutoUpdate', 'off', 'FontSize', 8)
                end
            end
            
            if row==1
                tit1 = sprintf('%s(%i)', data{sortInds(i), 'session'}{1}, data{sortInds(i), 'unit'});
                tit2 = data{sortInds(i), 'nucleus'}{1};
                title({tit1, tit2}, 'Interpreter', 'none', 'FontWeight', 'normal', 'FontSize', 8)
            end
        end
        
    else
        fprintf('WARNING! No spikes for %s unit %i\n', session, data{i, 'unit'})
    end
    
    % make new figure if necessary
    if mod(i, ncols)==0 || i==height(data)
        saveas(gcf, ['E:\lab_files\paper2\plots\rewards\reward_types\' num2str(figInd) '.png'])
        if i<height(data)
            figure('name', 'reward responses', 'color', 'white', 'position', figpos); hold on
            figInd = figInd + 1;
        end
    end
end












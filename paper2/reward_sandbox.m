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
tbl = table(init, init, init, nan(n,1), nan(n,1), nan(n,1), nan(n,1), 'VariableNames', ...
    {'response', 'predicted_full', 'predicted_noreward', 'dev_noreward', 'dev_full', 'mean', 'std'});

for i = 1:length(sessions)
    disp(i/length(sessions))
    % load ephys info
    neuralData = load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'neuralData', [sessions{i} '_neuralData.mat']), ...
        'unit_ids', 'spkRates', 'timeStamps');
    unit_ids = neuralData.unit_ids;
    load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'predictors', [sessions{i} '_predictors.mat']), 'predictors');
    rewardTimes = predictors{'reward_normal', 'data'}{1};
    X = rewardTimes + x;  % each row is the time grid for a single trial
    
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



















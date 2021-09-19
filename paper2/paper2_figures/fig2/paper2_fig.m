%% plot deviance explained for sequential models

data = getUnitInfo;
paper2_config;

%% example unit model predictions

%% show predictor traces surrounding reward delivery
close all;

% settings
sessions = {'200116_000', '191007_003', '201014_000'};
units = [30 3 193];

vars  = {'velocity', 'bodyAngle', 'paw1LH_x', 'paw2LF_x', 'paw3RF_x', 'paw4RH_x', 'jaw', 'whiskerAngle'};
names = {'velocity', 'body angle', 'left hind', 'left fore', 'right fore', 'right hind', 'jaw', 'whiskers'};
colors = [repmat(cfg.velColor, 6, 1); cfg.lickColor; cfg.wiskColor];

tlims = [-4 2];  % (s) time pre and post reward
offset = 4.5;  % (std) vertical offset for traces
figpos = [436.00 649.00 406.00 327.00];

vars = fliplr(vars); names = fliplr(names); colors = flipud(colors);

% sessions = data.session;
% units = data.unit;

for s = 1:length(sessions)
    try
        session = sessions{s};
        unit = units(s);

        % inits
        load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'predictors', [session '_predictors.mat']), 'predictors');
        load(['E:\lab_files\paper2\modelling\glms\sequential_glms\' session '_cell_' num2str(unit) '_glm.mat'], 'fitdata');
        validt = fitdata.t(~isnan(fitdata.yhat));
        validTLims = [min(validt) max(validt)];

    %     if mean(isnan(fitdata.yhat)) > 0
    %         disp(mean(isnan(fitdata.yhat)))
    %         load(['E:\lab_files\paper2\modelling\glms\sequential_glms\' session '_cell_' num2str(unit) '_glm.mat'], 'models');
    %         model = models{end, 'model'}{1};  % full model
    %         temp = exp(model.fit_preval) / dt;
    %         keyboard
    %     end

        figure('color', 'white', 'menubar', 'none', 'position', figpos);
        subplot(3,1,1:2); hold on
        rewardTimes = predictors{'reward_normal', 'data'}{1};
        rewardTimes = rewardTimes(rewardTimes>validTLims(1) & rewardTimes<validTLims(end));
        rewardTime = rewardTimes(fix(length(rewardTimes)/2));
        xlims = rewardTime + tlims;

        whiskerContactTimes = predictors{'whiskerContact', 'data'}{1};
        whiskerContactTimes = whiskerContactTimes(whiskerContactTimes>xlims(1) & whiskerContactTimes<xlims(2));
        t = predictors{vars{1}, 't'}{1};  % assumes all predictors have same time axis
        bins = t>=xlims(1) & t<=xlims(2);
        offsets = (0:(length(vars))) * offset;
        lickTimes = predictors{'lick', 'data'}{1};

        for i = 1:length(vars)
            sig = zscore(predictors{vars{i}, 'data'}{1}(bins)) + offsets(i);
            plot(t(bins), sig, 'LineWidth', 1.5, 'color', [colors(i, :) .6])
            text(xlims(1) - .01*diff(xlims), offsets(i), names{i}, 'HorizontalAlignment', 'right')

%             if strcmp(vars{i}, 'jaw')
%                 y = interp1(t(bins), sig, lickTimes);
%                 scatter(lickTimes, y, 20, cfg.lickColor, 'filled')
%             end
            if i==1; ymin = min(sig); end
        end

        % add line for events
        ylims = ylim; ylims(1) = ymin;
        ln = plot([rewardTime rewardTime], ylim, 'color', cfg.lickColor, 'LineWidth', 1);
        uistack(ln, 'bottom')
        text(rewardTime(1), ylims(2), 'reward', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')

        if ~isempty(whiskerContactTimes)
            ln = plot([whiskerContactTimes whiskerContactTimes]', ylim, 'color', cfg.wiskColor, 'LineWidth', 1);
            uistack(ln, 'bottom')
            text(whiskerContactTimes(1), ylims(2), 'whisker contact', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom')
        end

        set(gca, 'xlim', xlims, 'ylim', ylims, 'Visible', 'off')


        % firing rate vs. predicted firing rate
        subplot(3,1,3); hold on
        bins = fitdata.t>=xlims(1) & fitdata.t<=xlims(2);
        plot(fitdata.t(bins), fitdata.yRate(bins), 'color', [0 0 0 .4], 'LineWidth', 1.5)
        plot(fitdata.t(bins), fitdata.yhat(bins), 'color', cfg.predictionColor, 'LineWidth', 1.5)
        set(gca, 'XLim', xlims, 'Visible', 'off')
        legend('firing rate', 'predicted firing rate', 'autoupdate', 'off', 'FontSize', cfg.fontsize, 'box', 'off')

        % add time scale
        plot(xlims(1) + [0 0 .2], ylims(1) + [50 0 0], 'color', 'black', 'LineWidth', 1.5)
        yLims = ylim;
        text(xlims(1)+.1, yLims(1)-diff(yLims)*.05, '.2 second', ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'center')
        text(xlims(1)-diff(xlims)*.01, yLims(1)+25, '50 Hz', 'Rotation', 90, ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center')
        pause(.1);
        
        % save pngs to the option_figs folder
%         saveas(gcf, sprintf('E:\\lab_files\\paper2\\paper_figures\\option_figs\\prediction_unit_egs\\%s_unit_%i.png', session, unit))
        
        % save svgs for actual figgy
        saveas(gcf, ['E:\lab_files\paper2\paper_figures\matlab\model_example_' num2str(s) '.svg'])

    catch
        fprintf('ERROR WITH %s unit %i\n', session, unit)
    end
end





%% model performance for sequential models

% load sample glm
models = load(fullfile('E:\lab_files\paper2\modelling\glms\sequential_glms', ...
    sprintf('%s_cell_%i_glm.mat', data.session{1}, data.unit(1))));
modelNames = models.models.Properties.RowNames;
data.deviances = nan(height(data), length(modelNames));

% load all GLMs
for i = 1:height(data)
    fileName = fullfile('E:\lab_files\paper2\modelling\glms\sequential_glms', ...
    sprintf('%s_cell_%i_glm.mat', data.session{i}, data.unit(i)));
    if exist(fileName, 'file')
        models = load(fileName);
        data{i, 'deviances'} = models.models.dev';
    else
        fprintf('SKIPPING %s\n', fileName)
    end
end

%% plot
dmat = data.deviances';
close all; figure('color', 'white', 'position', [603.00 791.00 282.00 287.00], 'menubar', 'none')
bins = logical([0 1 1 1]);

dmatLog = dmat; dmatLog(dmatLog<=0) = nan; dmatLog = log(dmatLog);
barFancy(dmatLog(bins,:), 'showBars', false, 'showViolins', true, 'showErrorBars', false, ...
    'lineAtZero', false, 'levelNames', {modelNames(bins)'}, 'scatterColors', [0 0 0], ...
    'scatterSize', 10, 'ylabel', 'deviance explained', 'barWidth', .8, ...
    'colors', cfg.sequentialColors(bins,:), 'scatterCondColor', true, ...
    'comparisons', [1 2; 2 3], 'test', 'signrank', 'YLim', [-4.5 0]);
% barFancy(dmat, 'showBars', true, 'showViolins', false, 'showErrorBars', true, ...
%     'comparisons', [1 2; 2 3; 3 4], 'test', 'signrank', ...
%     'lineAtZero', true, 'levelNames', {modelNames'}, 'scatterColors', [0 0 0], ...
%     'scatterSize', 10, 'ylabel', 'deviance explained', 'barSeparation', 0, ...
%     'colors', cfg.sequentialColors, 'scatterCondColor', true, 'showScatter', false, ...
%     'barAlpha', 1, 'YLim', [0 .15], 'errorFunction', @(x) nanstd(x) / sqrt(length(x)));
limitticks
set(gca, cfg.axArgs{:});

% save
% set(gcf, 'Renderer', 'painters')  % ensures export is not rasterized
saveas(gcf, 'E:\lab_files\paper2\paper_figures\matlab\sequential_model_deviance.svg')
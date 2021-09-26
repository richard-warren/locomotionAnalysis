%% plot deviance explained for sequential models

data = getUnitInfo;
paper2_config;

%% load model performance

% load sample glm
models = load(fullfile('E:\lab_files\paper2\modelling\glms\sequential_glms', ...
    sprintf('%s_cell_%i_glm.mat', data.session{1}, data.unit(1))));
modelNames = models.models.Properties.RowNames;
data.deviances = nan(height(data), length(modelNames));
data.r2s = nan(height(data), length(modelNames));

% load all GLMs
for i = 1:height(data)
    fileName = fullfile('E:\lab_files\paper2\modelling\glms\sequential_glms', ...
        sprintf('%s_cell_%i_glm.mat', data.session{i}, data.unit(i)));
    if exist(fileName, 'file')
        models = load(fileName);
        data{i, 'deviances'} = models.models.dev';
        for j = 1:height(models.models)
            r2 = getModelR2(models.models{j, 'model'}{1}, models.fitdata);
            data.r2s(i, j) = r2;
        end
    else
        fprintf('SKIPPING %s\n', fileName)
    end
end

%% show predictor traces surrounding reward delivery
close all;

% settings
sessions = {'191007_003', '201014_000'};
units = [3 193];
predictorRows = 5;  % number of subplot rows occupied by predictors
exampleColors = lines(length(sessions));

modelsToShow = [2 3 4];
vars  = {'velocity', 'bodyAngle', 'paw1LH_x', 'paw2LF_x', 'paw3RF_x', 'paw4RH_x', '', 'whiskerAngle', 'jaw', 'reward_normal'};
names = {'velocity', 'body angle', 'left hind', 'left fore', 'right fore', 'right hind', '', 'whiskers', 'jaw', 'reward'};
colors = [repmat(cfg.velColor, 6, 1); .4 .4 .4; cfg.wiskColor; cfg.lickColor; cfg.lickColor];

tlims = [-4 2];  % (s) time pre and post reward
offset = 5;  % (std) vertical offset for traces
figpos = [436.00 507.00 406.00 469.00];


for s = 1:length(sessions)
    try
        session = sessions{s};
        unit = units(s);

        % inits
        load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'predictors', [session '_predictors.mat']), 'predictors');
        dmat = load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'designMatrices', 'sequential', [session '_designMatrix.mat']));
        load(['E:\lab_files\paper2\modelling\glms\sequential_glms\' session '_cell_' num2str(unit) '_glm.mat'], 'models', 'fitdata');
        validt = fitdata.t(~isnan(fitdata.yhat));
        validTLims = [min(validt) max(validt)];
        
        figure('color', 'white', 'menubar', 'none', 'position', figpos);
        subplot(predictorRows+length(modelsToShow), 1, 1:predictorRows); hold on
        rewardTimes = predictors{'reward_normal', 'data'}{1};
        rewardTimes = rewardTimes(rewardTimes>validTLims(1) & rewardTimes<validTLims(end));
        rewardTime = rewardTimes(fix(length(rewardTimes)/2));
        xlims = rewardTime + tlims;

        whiskerContactTimes = predictors{'whiskerContact', 'data'}{1};
        whiskerContactTimes = whiskerContactTimes(whiskerContactTimes>xlims(1) & whiskerContactTimes<xlims(2));
        t = predictors{vars{1}, 't'}{1};  % assumes all predictors have same time axis
        bins = t>=xlims(1) & t<=xlims(2);
        offsets = fliplr((0:(length(vars))) * offset);
        lickTimes = predictors{'lick', 'data'}{1};

        for i = 1:length(vars)
            if strcmp(vars{i}, '')
                scatter(repmat(mean(xlims), 1, 3), offsets(i) + [-offset 0 offset]*.3, 5, ...
                    colors(i,:), 'filled')
            elseif predictors{vars{i}, 'type'} == 'event'
                sig = 2*dmat.dmat{:, vars{i}}(bins,:) + offsets(i);
                sig(sig==offsets(i)) = nan;
                sig = [sig ones(size(sig,1),1)*offsets(i)];
                plot(t(bins), sig, 'LineWidth', 1.5, 'color', [colors(i, :) .6])
            else
                sig = zscore(predictors{vars{i}, 'data'}{1}(bins)) + offsets(i);
                plot(t(bins), sig, 'LineWidth', 1.5, 'color', [colors(i, :) .6])
            end
            text(xlims(1) - .01*diff(xlims), offsets(i), names{i}, 'HorizontalAlignment', 'right')
            if i==length(vars); ymin = min(sig(:)); end
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

        set(gca, 'xlim', xlims, 'YLim', ylims, 'Visible', 'off')
        
        % add title
        text(xlims(1), ylims(2), sprintf('EXAMPLE %i', s), 'Color', exampleColors(s,:), ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'center')


        % firing rate vs. predicted firing rate
        bins = fitdata.t>=xlims(1) & fitdata.t<=xlims(2);
        dt = t(2)-t(1);
        ymax = max(fitdata.yRate(bins));
        for j = 1:length(modelsToShow)
            subplot(predictorRows+length(modelsToShow), 1, predictorRows+j); hold on
            yhat = exp(models{modelsToShow(j), 'model'}{1}.fit_preval) / dt;
            plot(fitdata.t(bins), fitdata.yRate(bins), 'color', [0 0 0 .6], 'LineWidth', 1.5)
            plot(fitdata.t(bins), yhat(bins), 'color', cfg.sequentialColors(modelsToShow(j),:), 'LineWidth', 1.5)
            set(gca, 'XLim', xlims, 'Visible', 'off')
            text(xlims(1), ymax, cfg.sequentialModelNames{modelsToShow(j)}, 'HorizontalAlignment', 'left')
        end

        % add time scale
        legend('firing rate', 'predicted firing rate', ...
            'autoupdate', 'off', 'FontSize', cfg.fontsize, 'box', 'off')
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
dmat = data.deviances';
dmatLog = dmat; dmatLog(dmatLog<=0) = nan; dmatLog = log(dmatLog);
% close all
figure('color', 'white', 'position', [603.00 791.00 282.00 287.00], 'menubar', 'none')

% violin
% barFancy(dmatLog(bins,:), 'showBars', false, 'showViolins', true, 'showErrorBars', false, ...
%     'lineAtZero', false, 'levelNames', {modelNames(bins)'}, 'scatterColors', [0 0 0], ...
%     'scatterSize', 10, 'ylabel', 'deviance explained', 'barWidth', .8, ...
%     'colors', cfg.sequentialColors(bins,:), 'scatterCondColor', true, ...
%     'comparisons', [1 2; 2 3], 'test', 'signrank', 'YLim', [-4.5 0]);

% bar plot
% barFancy(dmat, 'showBars', true, 'showViolins', false, 'showErrorBars', true, ...
%     'comparisons', [1 2; 2 3; 3 4], 'test', 'signrank', ...
%     'lineAtZero', true, 'levelNames', {modelNames'}, 'scatterColors', [0 0 0], ...
%     'scatterSize', 10, 'ylabel', 'deviance explained', 'barSeparation', 0, ...
%     'colors', cfg.sequentialColors, 'scatterCondColor', true, 'showScatter', false, ...
%     'barAlpha', 1, 'YLim', [0 .15], 'errorFunction', @(x) nanstd(x) / sqrt(length(x)));

% box plot
newModelNames = cfg.sequentialModelNames;
colors = flipud(cfg.sequentialColors);
boxplot(dmat', newModelNames, ...
    'Symbol', '', 'Notch', 'on')
ylabel('fraction deviance explained')
h = findobj(gca,'Tag','Box');
for j = 1:length(h)
    patch(get(h(j),'XData'), get(h(j), 'YData'), colors(j,:), 'FaceAlpha',.5);
end
linesObs = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(linesObs, 'Color', [0 0 0], 'linewidth', 2);
lineObs = findobj(gcf, 'type', 'line', 'Tag', 'Lower Whisker'); set(lineObs, 'LineStyle', '-')
lineObs = findobj(gcf, 'type', 'line', 'Tag', 'Upper Whisker'); set(lineObs, 'LineStyle', '-')
set(gca, 'ylim', [-.05 .3], cfg.axArgs{:})
limitticks(true)
hold on

% add lines showing example units (need to run previous cell first!)
x = 1:length(newModelNames);
for i = 1:length(sessions)
    ind = find(strcmp(data.session, sessions{i}) & data.unit==units(i));
    y = dmat(:,ind);
    plot(x, y, '--', 'Color', [exampleColors(i,:) .6])
end

% save
saveas(gcf, 'E:\lab_files\paper2\paper_figures\matlab\sequential_model_deviance.svg')

%% scatter r2 against deviance explained

close all; figure('position', [504.00 724.00 343.00 285.00], 'color', 'white', 'menubar', 'none');
dev = data.deviances(:,4);
devLog = dev; devLog(devLog<0)=min(devLog(devLog>0))'; devLog = log(devLog);
r2 = data.r2s(:,4);
scatterHistoRick(dev, r2, 'colors', [0 0 0], 'scatAlpha', .1)
xlabel('deviance explained')
ylabel('R^2')

saveas(gcf, 'E:\lab_files\paper2\paper_figures\matlab\deviance_vs_r2.svg')

%% model performance on allen brain atlas

percentileLims = [10 90];

% inits
close all

dev = data.deviances(:, end);
% determine percentile based color map
clims = prctile(dev, percentileLims);
colors = interp1(clims, [cfg.sequentialColors(1,:); cfg.sequentialColors(end,:)], dev, 'nearest', 'extrap');

plotUnitsOnCcf(data, 'colors', colors)
saveas(gcf, 'E:\lab_files\paper2\paper_figures\matlab\ccf_deviance.svg')




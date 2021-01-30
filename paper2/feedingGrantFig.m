% tracking example frames for running and licking

session = '200709_000';
scatSize = 80;
smatScatSizes = [5 25];
contrast = [.1 .7];
trial = 2;
color = [0.9290 0.6940 0.1250];
pawColors = lines(4);
removeRows = [200 210];
runInd = 2875;
lickInd = 800;
featuresToShow = {'paw1LH_top', 'paw2LF_top', 'paw3RF_top', 'paw4RH_top', ...
                  'paw1LH_bot', 'paw2LF_bot', 'paw3RF_bot', 'paw4RH_bot', ...
                  'nose_top', 'nose_bot', 'tailMid_top', 'tailMid_bot', ...
                  'tailBase_top', 'tailBase_bot', 'jaw', 'ear', 'tongue'};

close all
for ind = [runInd lickInd]
    showSingleFrameTracking(session, 1, 'ind', ind, ...
        'contrastLims', contrast, 'addWiskCam', true, 'pawColors', pawColors, ...
        'otherColors', color, 'featuresToShow', featuresToShow, ...
        'removeRows', removeRows, 'mainSize', scatSize, 'trailingSizes', smatScatSizes);
    saveas(gcf, fullfile(getenv('OBSDATADIR'), 'data_transfer', sprintf('%stracking_frame_%i.png', session, ind)))
    
    % save face cam only
    vidWisk = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runWisk.mp4'));
    frameWisk = 255 - read(vidWisk, ind);
    imwrite(frameWisk, fullfile(getenv('OBSDATADIR'), 'data_transfer', sprintf('%sface_frame_%i.png', session, ind)))
end


%% show predictor traces surrounding reward delivery

session = '201016_000'; unit = 195;
rewardNum = 20;
vars = {'jaw',  'velocity', 'bodyAngle',  'paw4RH_x', 'whiskerAngle'};
names = {'jaw', 'velocity', 'body angle', 'paw',      'whiskers'};
tlims = [-6 4];  % (s) time pre and post reward
offset = 5;  % (std) vertical offset for traces
blue = [.5 .5 1];
figpos = [200.00 871.00 952.00 480.00];
glmColor = lines(1);

% inits
load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'predictors', [session '_predictors.mat']), 'predictors');
close all;

figure('color', 'white', 'menubar', 'none', 'position', figpos);
subplot(2,1,1); hold on
rewardTime = predictors{'reward_normal', 'data'}{1}(rewardNum);
xlims = rewardTime + tlims;
t = predictors{vars{1}, 't'}{1};  % assumes all predictors have same time axis
bins = t>=xlims(1) & t<=xlims(2);
offsets = (0:(length(vars))) * offset;
lickTimes = predictors{'lick', 'data'}{1};
vars = fliplr(vars); names = fliplr(names);

for i = 1:length(vars)
    sig = zscore(predictors{vars{i}, 'data'}{1}(bins)) + offsets(i);
    plot(t(bins), sig, 'LineWidth', 1.5, 'color', [glmColor .6])
    text(xlims(1) - .01*diff(xlims), offsets(i), names{i}, 'HorizontalAlignment', 'right')
    
    if strcmp(vars{i}, 'jaw')
        y = interp1(t(bins), sig, lickTimes);
        scatter(lickTimes, y, 40, blue, 'filled')
    end
end

% add line at reward time
yLims = ylim;
ln = plot([rewardTime rewardTime], ylim, 'color', blue, 'LineWidth', 2);
uistack(ln, 'bottom')

set(gca, 'xlim', xlims, 'ylim', ylim, 'Visible', 'off')


% firing rate vs. predicted firing rate
subplot(2,1,2); hold on
load(['E:\lab_files\paper2\modelling\glms\upper_lower_glms\' session '_cell_' num2str(unit) '_glm.mat'], 'fitdata');
bins = fitdata.t>=xlims(1) & fitdata.t<=xlims(2);

plot(fitdata.t(bins), fitdata.yRate(bins), 'color', [1 1 1]*.2)
plot(fitdata.t(bins), fitdata.yhat(bins), 'color', glmColor)
set(gca, 'XLim', xlims, 'Visible', 'off')
legend('firing rate', 'predicted firing rate')

saveas(gcf, 'Y:\loco\obstacleData\other\feedingGrantFigure\predictors', 'svg');


%% example cells, reward type comparison

% xcorr
close all

% settings
x = linspace(-2, 5, 200);  % x grid
units = struct();
units(1).session = '200624_000'; units(1).unit = 67;
units(2).session = '200722_000'; units(2).unit = 341;
colors = [0 0 0 1; 0 .6 0 .5; .8 0 0 .5];
corrDt = .01;
corrMaxLag = 1;  % (s)
kernelSz = .1;


figure('color', 'white', 'menubar', 'none', 'position', [200.00 660.00 341.00*length(units) 691.00]);
ncols = length(units);
maxLagSmps = corrMaxLag/corrDt;

for i = 1:ncols
   
    session = units(i).session; unit = units(i).unit;
    neuralData = load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'neuralData', [session '_neuralData.mat']));
        
    load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'predictors', [session '_predictors.mat']), 'predictors');
    rewardNormal = predictors{'reward_normal', 'data'}{1};
    rewardSurprise = predictors{'reward_surprise', 'data'}{1};
    rewardOmission = predictors{'reward_omission', 'data'}{1};

    conditions = {'normal', 'surprise', 'omission'};
    rewardTimes = {rewardNormal, rewardSurprise, rewardOmission};
    
    lickTimes = predictors{'lick', 'data'}{1};
    [lickRate, lickRateTimes] = getFiringRate(lickTimes, 'kernel', 'gauss', 'kernelSig', kernelSz, ...
        'fs', 100, 'tLims', [neuralData.timeStamps(1) neuralData.timeStamps(end)]);
    
    % compute firing rate with wider kernel
    unitInd = find(neuralData.unit_ids==unit);
    [fr, frTimes] = getFiringRate(neuralData.spkTimes{unitInd}, 'kernel', 'gauss', 'kernelSig', kernelSz);
    
    % times all times and firing rate into table
    unitData = table({lickRateTimes, frTimes}', ...
                     {lickRate, fr}', ...
                     'VariableNames', {'t', 'data'}, ...
                     'RowNames', {'lick rate', 'firing rate'});
    
    for row = 1:height(unitData)
        subplot(height(unitData)+1, ncols, i+(row-1)*ncols); hold on

        for j = 1:3  % for each of 3 response types
            rewardBins = rewardTimes{j}>=frTimes(1) & rewardTimes{j}<=frTimes(end);  % bins where unit is not nan
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
        set(gca, 'ylim', yLims, 'xlim', [x(1) x(end)])

        if i==1
            ylabel(unitData.Properties.RowNames{row});
            if row==height(unitData); xlabel('time after reward (s)'); end
            if row==1
                temp = [];
                for k = 1:3; temp(k) = plot([nan nan], 'color', colors(k,:), 'LineWidth', 2); end % create dummy lines
                legend(temp, conditions, 'Location', 'best', 'Box', 'off', ...
                    'AutoUpdate', 'off', 'FontSize', 8)
            end
        end
    end
    
    % xcorr
    subplot(height(unitData)+1, ncols, i+2*ncols); hold on
    t = neuralData.timeStamps(~isnan(neuralData.spkRates(unitInd,:)));
    t = t(1) : corrDt : t(end);  % interpole onto common .01 s time grid
    lickRateInterp = zscore(interp1(lickRateTimes, lickRate, t));
    frInterp = zscore(interp1(neuralData.timeStamps, neuralData.spkRates(unitInd,:), t));
    [corrLick, lags] = xcorr(frInterp, lickRateInterp, maxLagSmps, 'unbiased');
    plot(lags*corrDt, corrLick, 'LineWidth', 2, 'color', 'black')
    ylims = ylim;
    plot([0 0], ylims, 'color', [1 1 1]*.4);
    [~, maxInd] = max(abs(corrLick));
    scatter(lags(maxInd)*corrDt, corrLick(maxInd), 50, 'black', 'filled');
    xlabel('<- brain leads (lag [s]) body leads ->');
    ylabel('correlation');
end

saveas(gcf, 'Y:\loco\obstacleData\other\feedingGrantFigure\reward_types', 'svg');


%% example cells, lick rate comparison
% lick rate, broken down by integral
% firing rate, broken down by integral
% lick rate, by reward condition
% firing rate, by reward condition
% xcorr
close all

% settings
units = struct();
units(1).session = '200624_000'; units(1).unit = 67;
units(2).session = '200722_000'; units(2).unit = 341;

% fastigialBins = strcmp(data.nucleus, 'fastigial');
% units = struct('session', data.session(fastigialBins), ...
%                'unit', num2cell(data.unit(fastigialBins)));
% units = units(61:end);

x = linspace(-2, 5, 200);  % x grid
nbins = 3;  % number of lick response bins
colors = [lines(1)*.5; lines(1)];
corrDt = .01;
corrMaxLag = 1;  % (s)
kernelSz = .1;
highLowPercentiles = [20 80];
conditions = {'low lick rate', 'high lick rate'};


figure('color', 'white', 'menubar', 'none', 'position', [200.00 660.00 341.00*length(units) 691.00]);
ncols = length(units);
maxLagSmps = corrMaxLag/corrDt;

for i = 1:ncols
   
    session = units(i).session; unit = units(i).unit;
    neuralData = load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'neuralData', [session '_neuralData.mat']));
    unitInd = find(neuralData.unit_ids==unit);
    [fr, frTimes] = getFiringRate(neuralData.spkTimes{unitInd}, 'kernel', 'gauss', 'kernelSig', kernelSz);
    
    load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'predictors', [session '_predictors.mat']), 'predictors');
    rewardTimes = predictors{'reward_normal', 'data'}{1};
    rewardTimes = rewardTimes(rewardTimes>frTimes(1) & rewardTimes<frTimes(end));
    
    lickTimes = predictors{'lick', 'data'}{1};
    [lickRate, lickRateTimes] = getFiringRate(lickTimes, 'kernel', 'gauss', 'kernelSig', kernelSz, ...
        'fs', 100, 'tLims', [neuralData.timeStamps(1) neuralData.timeStamps(end)]);
    lickResp = interp1(lickRateTimes, lickRate, rewardTimes + x);
    
    % determine high and low lick rate trials
    lickAmp = mean(lickResp,2);
    binEdges = prctile(lickAmp, [0 highLowPercentiles 100]);
    [~, ~, binNums] = histcounts(lickAmp, binEdges);
    binNums(binNums==2) = nan;
    binNums(binNums==3) = 2;
    
    % times all times and firing rate into table
    unitData = table({lickRateTimes, frTimes}', ...
                     {lickRate, fr}', ...
                     'VariableNames', {'t', 'data'}, ...
                     'RowNames', {'lick rate', 'firing rate'});
    
    for row = 1:height(unitData)
        subplot(height(unitData)+1, ncols, i+(row-1)*ncols); hold on

        for j = 1:2  % for each of 2 lick conditions
            rewardBins = binNums==j;  % bins where unit is not nan
            resp = interp1(unitData{row, 't'}{1}, ...
                           unitData{row, 'data'}{1}, ...
                           rewardTimes(rewardBins) + x);
            mn = nanmean(resp,1);
            dev = nanstd(resp,1);
            patch([x fliplr(x)], [mn+dev fliplr(mn-dev)], colors(j,:), ...
                'EdgeColor', 'none', 'FaceAlpha', .1)
            plot(x, mn, 'Color', colors(j,:), 'LineWidth', 3)
        end
        yLims = ylim;
        plot([0 0], yLims, 'Color', [0 0 0 .5])  % line at x=0
        set(gca, 'ylim', yLims, 'xlim', [x(1) x(end)])

        if i==1
            ylabel(unitData.Properties.RowNames{row});
            if row==height(unitData); xlabel('time after reward (s)'); end
            if row==1
                temp = [];
                for k = 1:2; temp(k) = plot([nan nan], 'color', colors(k,:), 'LineWidth', 2); end % create dummy lines
                legend(temp, conditions, 'Location', 'northeast', 'Box', 'off', ...
                    'AutoUpdate', 'off', 'FontSize', 8)
            end
        end
        if row==1; title(sprintf('%s (%i)', session, unit), 'Interpreter', 'none'); end
    end
    
    % xcorr
    subplot(height(unitData)+1, ncols, i+2*ncols); hold on
    t = neuralData.timeStamps(~isnan(neuralData.spkRates(unitInd,:)));
    t = t(1) : corrDt : t(end);  % interpole onto common .01 s time grid
    lickRateInterp = zscore(interp1(lickRateTimes, lickRate, t));
    frInterp = zscore(interp1(neuralData.timeStamps, neuralData.spkRates(unitInd,:), t));
    [corrLick, lags] = xcorr(frInterp, lickRateInterp, maxLagSmps, 'unbiased');
    plot(lags*corrDt, corrLick, 'LineWidth', 2, 'color', 'black')
    ylims = ylim;
    plot([0 0], ylims, 'color', [1 1 1]*.4);
    [~, maxInd] = max(abs(corrLick));
    scatter(lags(maxInd)*corrDt, corrLick(maxInd), 50, 'black', 'filled');
    xlabel('<- brain leads (lag [s]) body leads ->');
    ylabel('correlation');
end

saveas(gcf, 'Y:\loco\obstacleData\other\feedingGrantFigure\high_low_licks', 'svg');


%% histo panels

units = struct();
units(1).session = '200624_000'; units(1).unit = 67;
units(2).session = '200722_000'; units(2).unit = 341;


ccf = loadCCF();

%%
close all
figure('color', 'white', 'menubar', 'none', 'position', [757.00 583.00 313.00 646.00]);
colors = repelem(lines(3), 2, 1);
args = {'apGrid', ccf.ap, 'dvGrid', ccf.dv, 'mlGrid', ccf.ml};
planeColor = lines(1);

subplot(2,1,1); hold on
plotLabels3D(ccf.coarseLabels, 'downSampling', 7, 'colors', zeros(3,3), ...
    'surfArgs', {'FaceAlpha', .04, 'EdgeColor', 'none'}, args{:})
plotLabels3D(ccf.labels, 'colors', colors, 'downSampling', 2, ...
    'surfArgs', {'EdgeColor', 'none', 'FaceAlpha', .8}, args{:})
scatter3(data.ccfMm(:,1), data.ccfMm(:,2), data.ccfMm(:,3), ...
    4, 'black', 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeColor', 'none');
set(gca, 'XTick', [], 'YTick', [], 'ZTick', [])
view(45, 30)

% draw plan at constant coronal location
xlims = xlim; zlims = zlim;
ap = 10;
patch([xlims(1) xlims(2) xlims(2) xlims(1)], ...
      [ap ap ap ap],...
      [zlims(1) zlims(1) zlims(2) zlims(2)], ...
      planeColor, 'EdgeColor', planeColor, 'FaceAlpha', .1, 'linewidth', 3)  % (ml, ap, dv)

subplot(2,1,2); hold on
plotLabels2D(ccf.labels, 'dim', 'ap', 'colors', colors, args{:})
set(gca, 'visible', 'off')
xlims = xlim; ylims = ylim;
patch([xlims(1) xlims(2) xlims(2) xlims(1)], ...
      [ylims(1) ylims(1) ylims(2) ylims(2)], ...
      planeColor, 'EdgeColor', planeColor, 'FaceColor', 'none', 'linewidth', 1)  % (ml, ap, dv)


% set(gca, 'ylim', [2.5 5])
scatter(data.ccfMm(:,1), data.ccfMm(:,3), ...
        10, 'black', 'filled', 'MarkerFaceAlpha', .6, 'MarkerEdgeColor', 'none');

set(gcf, 'Renderer', 'painters')
saveas(gcf, 'Y:\loco\obstacleData\other\feedingGrantFigure\histo', 'svg');







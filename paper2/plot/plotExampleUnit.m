function plotExampleUnit(session, unit)
% plot unit firing rate and predictor traces around reward, and a handful
% of PSTHs to show tuning to different predictors

% settings
paper2_config;  % creates 'cfg' variable
s.predictors = {'velocity', 'paw4RH_x', 'jaw'};
s.predictorNames = {'velocity', 'right fore', 'jaw'};
s.psths = {'velocity', 'reward_normal', 'whiskerContact', 'paw4RH_stride'};
s.psthNames = {'velocity', 'reward', 'whisker contact', 'phase'};
s.xlabels = {'velocity (m/s)', 'time from reward (s)', 'time from whisker contact (s)', 'fraction of stride'};

s.tlims = [-4 2];  % (s) time pre and post reward
s.offset = 5;  % (std) vertical offset for traces

s.predictorColors = [cfg.velColor; cfg.velColor; cfg.lickColor];
s.psthColors = [cfg.velColor; cfg.lickColor; cfg.wiskColor; cfg.velColor];
s.psthxlims = [nan nan; -.4 1; -.2 .4; 0 1];  % nan for auto determination
s.nscatters = 1000;  % scatter points to include for continuous variables






% inits
load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'predictors', [session '_predictors.mat']), 'predictors');
fr = load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'neuralData', [session '_neuralData.mat']), ...
    'unit_ids', 'spkRates', 'timeStamps');
unitInd = find(fr.unit_ids==unit);
fr.fr = fr.spkRates(unitInd, :);
load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'responses', [session '_responses.mat']), 'responses');

figure('color', 'white', 'menubar', 'none', 'position', [7.00 1156.00 1265.00 139.00])
ncols = length(s.psths)+2;
subplot(1, ncols, 1:2); hold on

validt = fr.timeStamps(~isnan(fr.fr));
validTLims = [min(validt) max(validt)];


rewardTimes = predictors{'reward_normal', 'data'}{1};
rewardTimes = rewardTimes(rewardTimes>validTLims(1) & rewardTimes<validTLims(end));
rewardTime = rewardTimes(fix(length(rewardTimes)/2));  % pick middle reward
xlims = rewardTime + s.tlims;

whiskerContactTimes = predictors{'whiskerContact', 'data'}{1};
whiskerContactTimes = whiskerContactTimes(whiskerContactTimes>xlims(1) & whiskerContactTimes<xlims(2));
t = predictors{s.predictors{1}, 't'}{1};  % assumes all predictors have same time axis
bins = t>=xlims(1) & t<=xlims(2);
offsets = fliplr((0:length(s.predictors)) * s.offset);

% predictor traces
for i = 1:length(s.predictors)
    sig = zscore(predictors{s.predictors{i}, 'data'}{1}(bins)) + offsets(i);
    plot(t(bins), sig, 'LineWidth', 1.5, 'color', [s.predictorColors(i, :) .6])
    text(xlims(1) - .01*diff(xlims), offsets(i), s.predictorNames{i}, 'HorizontalAlignment', 'right')
end

% firing rate
bins = fr.timeStamps>=xlims(1) & fr.timeStamps<=xlims(2);
sig = zscore(fr.fr(bins));
ymin = min(sig);
plot(fr.timeStamps(bins), sig, 'LineWidth', 1.5, 'color', [.15 .15 .15])
text(xlims(1) - .01*diff(xlims), 0, 'firing rate', 'HorizontalAlignment', 'right')

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


% add time scale
% TODO: add firing rate
plot(xlims(1) + [0 .2], ylims(1) + [0 0], 'color', 'black', 'LineWidth', 1.5)
yLims = ylim;
text(xlims(1)+.1, yLims(1)-diff(yLims)*.05, '.2 second', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'center')

for i = 1:length(s.psths)
    subplot(1, ncols, 2+i); hold on
    title(s.psthNames{i})
    
    type = responses{s.psths{i}, 'type'};
    resp = responses{s.psths{i}, 'response'}{1};
    xlims = responses{s.psths{i}, 'xLims'};
    
    if type == 'continuous'
        x = linspace(xlims(1), xlims(2), size(resp, 1));
        
        % scatters
        frInterp = interp1(fr.timeStamps, fr.fr, t);
        validbins = ~isnan(frInterp);
        frInterp = frInterp(validbins);
        temp = predictors{s.psths{i}, 'data'}{1}(validbins);
        inds = randperm(length(frInterp), s.nscatters);
        scatter(temp(inds), frInterp(inds), 5, s.psthColors(i,:), 'filled', 'markerfacealpha', .25)
        
        % density
        density = responses{s.psths{i}, 'density'}{1};
        ymax = max(frInterp(inds));
        density = density * (ymax / max(density));
        dens = fill([xlims(1) x xlims(2) xlims(1)]', ...
            [0 density 0 0]', [0 0 0], ...
            'EdgeColor', [1 1 1]*.6, 'FaceAlpha', .1);
        uistack(dens, 'bottom')
        
        % moving average
        plot(x, resp(:, unitInd), 'LineWidth', 3, 'color', s.psthColors(i,:))
    else
        resp = resp(:, :, unitInd);
        mn = nanmean(resp, 1);
        st = nanstd(resp, 1);
        x = linspace(xlims(1), xlims(2), size(resp, 2));
        
        patch([x fliplr(x)], [(-st+mn) fliplr(st+mn)], s.psthColors(i,:), ...
            'FaceAlpha', .25, 'EdgeColor', 'none')           % shaded error bars
        plot(x, mn, 'color', s.psthColors(i,:), 'lineWidth', 2);  % mean
        
        if type == 'event'
            plot([0 0], ylim, 'color', [.4 .4 .4], 'LineWidth', 1, 'color', [0 0 0 .4]);
        end
    end
    
    if ~any(isnan(s.psthxlims(i,:))); xlims = s.psthxlims(i,:); end
    set(gca, 'XLim', xlims, cfg.axArgs{:})
    xlabel(s.xlabels{i})
    
    if i==1; ylabel('firing rate'); end
end






















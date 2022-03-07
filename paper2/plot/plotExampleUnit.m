function plotExampleUnit(session, unit, groups, groupDevExplained, varargin)
% plot unit firing rate and predictor traces around reward, and a handful
% of PSTHs to show tuning to different predictors // `groupDevExplained` is
% a 2 X ngroups matrix of the lower and upper bounds of deviance explained
% for each group in `groups`

% settings
paper2_config;  % creates 'cfg' variable
s.predictors = {'velocity', 'paw4RH_x', 'jaw'};
s.predictorNames = {'velocity', 'right fore', 'jaw'};

s.psths = {'velocity', 'paw4RH_stride', 'reward_normal', 'lick', 'whiskerContact'};
s.psthNames = {'velocity', 'phase', 'reward', 'lick', 'whisker contact'};
s.xlabels = {'velocity (m/s)', 'fraction of stride', 'time from reward (s)', 'time from lick (s)', 'time from whisker contact (s)'};
s.psthxlims = [nan nan; 0 1; -.4 1; -.2 .2; -.2 .4];  % nan for auto determination
s.psthColors = [cfg.velColor; cfg.velColor; cfg.lickColor; cfg.lickColor; cfg.wiskColor];

s.predictorColors = [cfg.velColor; cfg.velColor; cfg.lickColor];
s.tlims = [-4 2];  % (s) time pre and post reward
s.offset = 5;  % (std) vertical offset for traces
s.nscatters = 1000;  % scatter points to include for continuous variables
s.showDeviance = true;

% a hack to stack subplots on top of eachother in an other loop
s.subplotRows = 1;
s.subplotColIdx = 1; 
s.hideTraceText = false;
s.hideTitles = false;
s.hideX = false;

s.title = '';




% inits
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'predictors', [session '_predictors.mat']), 'predictors');
fr = load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'neuralData', [session '_neuralData.mat']), ...
    'unit_ids', 'spkRates', 'timeStamps');
unitInd = find(fr.unit_ids==unit);
fr.fr = fr.spkRates(unitInd, :);
load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'responses', [session '_responses.mat']), 'responses');

if s.subplotRows == 1; figure('color', 'white', 'menubar', 'none', 'position', [7.00 1156.00 1265.00 139.00]); end
ncols = length(s.psths) + 2 + s.showDeviance;


% DEVIANCE EXPLAINED
if s.showDeviance
    subplot(s.subplotRows, ncols, 3 + ncols*(s.subplotColIdx-1)); hold on
    x = 1:length(groups);
    plot([x; x], groupDevExplained, 'color', [.15 .15 .15])
    scatter(x, groupDevExplained(1,:), 10, cfg.upperLowerColors(1,:), 'filled');
    scatter(x, groupDevExplained(2,:), 10, cfg.upperLowerColors(2,:), 'filled');
    set(gca, 'xtick', x, 'XTickLabel', groups, 'XTickLabelRotation', 30, ...
        cfg.axArgs{:}, 'FontSize', 6)
    limitticks(true)
end

% SAMPLE TIME TRACES
subplot(s.subplotRows, ncols, (1:2) + ncols*(s.subplotColIdx-1)); hold on

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
    if ~s.hideTraceText
        text(xlims(1) - .01*diff(xlims), offsets(i), s.predictorNames{i}, 'HorizontalAlignment', 'right')
    end
end

% firing rate
bins = fr.timeStamps>=xlims(1) & fr.timeStamps<=xlims(2);
sig = zscore(fr.fr(bins));
ymin = min(sig);
plot(fr.timeStamps(bins), sig, 'LineWidth', 1.5, 'color', [.15 .15 .15])
if ~s.hideTraceText; text(xlims(1) - .01*diff(xlims), 0, 'firing rate', 'HorizontalAlignment', 'right'); end

% add line for events
ylims = ylim; ylims(1) = ymin;
ln = plot([rewardTime rewardTime], ylim, 'color', cfg.lickColor, 'LineWidth', 1);
uistack(ln, 'bottom')
if ~s.hideTraceText; text(rewardTime(1), ylims(2), 'reward', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom'); end

if ~isempty(whiskerContactTimes)
    ln = plot([whiskerContactTimes whiskerContactTimes]', ylim, 'color', cfg.wiskColor, 'LineWidth', 1);
    uistack(ln, 'bottom')
    if ~s.hideTraceText
        text(whiskerContactTimes(1), ylims(2), 'whisker contact', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom')
    end
end

set(gca, 'xlim', xlims, 'YLim', ylims, 'Visible', 'off')

% add title
if ~isempty(s.title)
    text(mean(xlims), ylims(1), s.title, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'FontWeight', 'bold', 'Interpreter', 'none')
end


% PSTHS

% add time scale
% TODO: add firing rate
plot(xlims(1) + [0 .5], ylims(1) + [0 0], 'color', 'black', 'LineWidth', 1.5)
yLims = ylim;
if ~s.hideTraceText
    text(xlims(1)+.1, yLims(1)-diff(yLims)*.05, '.5 second', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'center')
end

for i = 1:length(s.psths)
    subplot(s.subplotRows, ncols, 2+s.showDeviance+i + ncols*(s.subplotColIdx-1)); hold on
    if ~s.hideTitles; title(s.psthNames{i}); end
    
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
            ylims = ylim;
            plot([0 0], ylim, 'color', [.4 .4 .4], 'LineWidth', 1, 'color', [0 0 0 .4]);
            set(gca, 'ylim', ylims);
        end
    end
    
    if ~any(isnan(s.psthxlims(i,:))); xlims = s.psthxlims(i,:); end
    set(gca, 'XLim', xlims, cfg.axArgs{:})
    if ~s.hideX; xlabel(s.xlabels{i}); end
    limitticks
    if s.hideX; set(gca, 'XTickLabel', []); end
    
    if i==1; ylabel('firing rate'); end
end






















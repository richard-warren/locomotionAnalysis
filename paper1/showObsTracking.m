function showObsTracking(session, varargin)

% show wheel vs. obstacle velocity for an example session to demonstrate
% how obstacle successfully tracks wheel velocity

% temp

% settings
s.figPos = [662.00 425.00 362.00 305.00];
s.numTrials = 10;
s.velTime = .05;
s.wheelColor = [.2 .2 .2];
s.obsColor = [0.7373 0.4902 0.7098];
s.waterColor = 'blue';
s.offset = .6;  % (m/s)
s.xScale = 1;  % (s) x scale bar
s.yScale = .8;  % (m/s) y scale bar
s.minTime = 3;  % how many seconds to trim off the left side of plot (when mouse is licking water from previous reward)
s.durationPercentile = 50; % only take trials with duration less that durationPercentile
s.smpsPerTrial = 1000;  % resample so each trial only has this many samples // this allows the vector to be preserved and keeps file sizes low


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'wheelPositions', 'wheelTimes', 'obsPositions', 'obsTimes', ...
    'obsOnTimes', 'obsOffTimes', 'rewardTimes')
wheelVel = getVelocity(wheelPositions, s.velTime, 1/median(diff(wheelTimes)));
obsVel = getVelocity(obsPositions, s.velTime, 1/median(diff(obsTimes)));
figure('name', 'obsTracking', 'Color', 'white', 'MenuBar', 'none', ...
    'Position', s.figPos, 'inverthardcopy', 'off'); hold on

% determine when obstacle is on
onInds = knnsearch(obsTimes', obsOnTimes);
offInds = knnsearch(obsTimes', obsOffTimes);
isObsOn = zeros(size(obsTimes));
isObsOn(onInds) = true;
isObsOn(offInds) = -1;
isObsOn = logical(cumsum(isObsOn));

% select trials
durations = diff(rewardTimes);
maxDuration = prctile(durations, s.durationPercentile);
validTrials = find(durations<maxDuration);
trials = validTrials(1:s.numTrials);




for i = 1:s.numTrials
    
    % get wheel data
    bins = wheelTimes>=rewardTimes(trials(i)) & wheelTimes<=rewardTimes(trials(i)+1);
    w = wheelVel(bins);
    wt = wheelTimes(bins);
    
    % get obstacle data
    bins = obsTimes>=rewardTimes(trials(i)) & obsTimes<=rewardTimes(trials(i)+1);
    o = obsVel(bins);
    ot = obsTimes(bins);
    o(~isObsOn(bins)) = nan;  % mask regions where the obstacle is not on
    
    % add vertical offset
    w = w + (i-1)*s.offset;
    o = o + (i-1)*s.offset;
    
    % interpolate over smaller grid
    wt_interp = linspace(wt(1), wt(end), s.smpsPerTrial);
    w = interp1(wt, w, wt_interp);
    ot_interp = linspace(ot(1), ot(end), s.smpsPerTrial);
    o = interp1(ot, o, ot_interp);
    
    
    plot(wt_interp-wt_interp(1), w, 'LineWidth', 1, 'Color', s.wheelColor); hold on
    plot(ot_interp-ot_interp(1), o, 'LineWidth', 2, 'Color', s.obsColor)
    
    scatter(range(wt_interp), w(end), 20, s.waterColor, 'filled')  % add scatter at point of water reward
end

xLims = [s.minTime, max(durations(trials))];
set(gca, 'xlim', xLims, 'visible', 'off', 'ylim', [-.1, max(w)])

yMin = min(get(gca, 'ylim'));
x = [xLims(2)-s.xScale, xLims(2)];
y = [yMin, yMin + s.yScale];
plot([x(1) x(2) x(2)], [y(1) y(1) y(2)], 'lineWidth', 3, 'color', get(gca, 'xcolor'))  % add scale bar
text(mean(x), y(1), sprintf('%i seconds', s.xScale), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top')
text(x(2), y(2), sprintf('%.1f m/s', s.yScale), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')


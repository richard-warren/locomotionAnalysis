function plotPSTH(spkTimes, events, varargin)

% plot PSTH with instantaneous firing rate on the top and raster on the
% bottom // give spike times and events vector // events can also be nX2
% matrix of 'epoch' start and stop times, in which case epochs are 
% stretched over common x axis

% todo: restrict bin durations?



% settings
% --------

% x axis
s.eventLims = [-.25 .5];  % (s) x limits for events
s.epochLims = [-.25 1.25];  % (fraction of epoch) x limits for epochs
s.binNum = 500;  % number of points on the axis

% raster
s.scatSize = 3.5;

% instantaneous firing rate kernel
s.kernelRise = .005;     % (s) rise for double exponential kernel
s.kernelFall = .02;      % (s) fall for double exponential kernel
s.kernelSig = .02;       % (s) if a gaussian kernel is used
s.kernel = 'doubleExp';  % 'gauss', or 'doubleExp'

% condition info
s.conditions = ones(size(events,1),1);  % vector of trial condition numbers
s.color = [];  % colors for each condition
s.conditionNames = {};

% other
s.epochDurationLims = [];  % (percentiles) only include epochs of duration within these limits
s.xlabel = '';
s.removeNoSpikeTrials = false;  % whether to automatically remove trials with no spikes
s.maxEpochs = 2000;  % if using epoch events and there are more than maxEpochs, limit to epochs of central duration // set to 0 to bypass


% initializations
% ---------------
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end % reassign settings passed in varargin
conditions = unique(s.conditions);

if isempty(s.color)
    if length(conditions)>1
        s.color = lines(length(conditions));
    else
        s.color = [1 1 1]*.15;  % gray if only one condition
    end
end

isEpoch = size(events,2)==2;  % whether event or epoch
[spkRates, spkRatesTimes] = getFiringRate(spkTimes, 'fs', 2000, ...
    'kernel', s.kernel, 'kernelRise', s.kernelRise, 'kernelFall', s.kernelFall, 'kernelSig', s.kernelSig);
if isEpoch; xLims = s.epochLims; else; xLims = s.eventLims; end
x = linspace(xLims(1), xLims(2), s.binNum);

% if too many epochs find the elements of central duration
if isEpoch && s.maxEpochs
    durations = diff(events,1,2);
    if size(events,1)>s.maxEpochs
        middleIndsSorted = floor(length(durations)/2 - s.maxEpochs*.5) + (1:s.maxEpochs);
        [~, sortInds] = sort(durations);
        inds = sortInds(middleIndsSorted)';  % !!! should maybe pick random intead of central-duration epochs
        
        events = events(inds,:);
        s.conditions = s.conditions(inds);
    end
end


% get responses
% -------------
responseSpks = cell(length(events),1);  % vector of spike times for every trial
if ~isEpoch
    responseRate = interp1(spkRatesTimes, spkRates, x+events);  % compute inst firing rate for all trials in one shot
else
    responseRate = nan(length(events), s.binNum);
end

for i = 1:length(events)
    
    % events
    if ~isEpoch
        responseSpks{i} = spkTimes(spkTimes>=(events(i)+xLims(1)) & spkTimes<=(events(i)+xLims(2))) - events(i);
        
    % epochs
    else
        epoch = range(events(i,:))*s.epochLims + events(i,1);  % start and end time for epoch
        bins = spkRatesTimes>=epoch(1) & spkRatesTimes<=epoch(2);
        if any(bins)
            responseRate(i,:) = interp1(spkRatesTimes(bins), spkRates(bins), ...
                linspace(epoch(1), epoch(2), s.binNum), 'linear');
            responseRate(i,:) = fillmissing(responseRate(i,:), 'linear', 'EndValues', 'nearest');
            
            trialSpks = spkTimes(spkTimes>=epoch(1) & spkTimes<=epoch(2));
            responseSpks{i} = (trialSpks-events(i,1)) * (range(xLims) / range(epoch));  % rescale and shift onto new x axis
        end 
    end
end

% remove trials with no spikes
if s.removeNoSpikeTrials
    notEmpty = ~cellfun(@isempty, responseSpks);
    
    events = events(notEmpty,:);
    responseSpks = responseSpks(notEmpty);
    responseRate = responseRate(notEmpty,:);
    s.conditions = s.conditions(notEmpty);
end

% remove outlier epochs
if isEpoch && ~isempty(s.epochDurationLims)
    durations = diff(events,1,2);
    durationLims = prctile(durations, s.epochDurationLims);
    validDurations = durations>=durationLims(1) & durations<= durationLims(2);
    
    events = events(validDurations,:);
    responseSpks = responseSpks(validDurations);
    responseRate = responseRate(validDurations,:);
    s.conditions = s.conditions(validDurations);
end



% plot
% ----
figure('color', 'white', 'menubar', 'none', 'position', [1380.00 157.00 331.00 847.00])
meanLines = nan(length(conditions), 1);

% firing rate
subplot(3,1,1); hold on
for i = 1:length(conditions)
    bins = s.conditions==conditions(i);
    condMean = nanmean(responseRate(bins,:),1);
    condStd = nanstd(responseRate(bins,:),1);
    
    meanLines(i) = plot(x, condMean, 'color', s.color(i,:), 'LineWidth', 3);
    patch([x(1) x fliplr(x)], ...
        [condMean(1)-condStd(1), condMean+condStd fliplr(condMean-condStd)], ...
        s.color(i,:), 'facealpha', .2, 'edgecolor', 'none')
end

yLims = max(get(gca, 'YLim'), 0);
addVerticalLines(isEpoch, xLims, yLims);
set(gca, 'xlim', xLims, 'ylim', yLims)
if ~isempty(s.conditionNames)
    legend(meanLines, s.conditionNames, 'Location', 'best', 'Box', 'off')
end
ylabel('firing rate (hz)')

% histogram
subplot(3,1,2:3); hold on
[conditionsSorted, sortInds] = sort(s.conditions);
responseSpksSorted = responseSpks(sortInds);
xScat = cat(1, responseSpksSorted{:});
yScat = repelem(1:size(events,1), cellfun(@length, responseSpks));
spkConditions = repelem(conditionsSorted, cellfun(@length, responseSpks));
scatter(xScat, yScat, s.scatSize, s.color(spkConditions,:), 'filled', ...
    'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', .6)

yLims = [1 size(events,1)];
addVerticalLines(isEpoch, xLims, yLims);
set(gca, 'xlim', xLims, 'ylim', yLims, 'YTick', [], 'TickDir', 'out', 'YDir', 'normal')
ylabel('spikes')
if ~isempty(s.xlabel); xlabel(s.xlabel, 'Interpreter', 'none'); end



% functions
% ---------

function addVerticalLines(isEpoch, xLims, yLims)
    if isEpoch
        if xLims(1)<0; plot([0 0], yLims, 'LineWidth', 1, 'Color', [0 0 0 .4]); end
        if xLims(2)>1; plot([1 1], yLims, 'LineWidth', 1, 'Color', [0 0 0 .4]); end
    else
        plot([0 0], yLims, 'LineWidth', 1, 'Color', [0 0 0 .4]);
    end
end

end


function plotPSTH2(session, cellNum, eventTimes, opts)

% to do: interp all cells within a session at the same time (to speed
% things up) // plot multiple thing on top of eachother, which would let me
% get rid of bootstrappinfg too... // add vertical line option for epoch
% mode, eg to show swing and stance // accept previously computed data, and
% replot only, mwahahahaha // number plots // add buffer to begining and
% end of epoch plots

% !!! need to document


% settings
s.normalize = false;  % whether to z score data
s.yLimsNormalized = [-2.5 2.5]; % plot mean +- yLimsNormalized*std
s.xLims = [-.1 .5]; % limits for x axis
s.plotTrials = 0; % if >0, plot individual trial firing rates in the background
s.xGridLength = 100; % number of points per trial for epoch interpolation
s.showSampleSize = true;  % whether to add number of trials in parentheses to the plot
s.conditionNames = {};
s.conditionColors = [];
s.errorFcn = @(x) nanstd(x); % function for error bars // set to false if you don't want any error bars!
% s.errorFcn = @(x) nanstd(x)/sqrt(size(x,1)); % function for error bars
s.colors = 'hsv'; % color scheme // can be specified either as a string, or as an nX3 color map, where n is number of conditions
s.trialAlpha = .1;  % transparency of individual trial lines
s.errorAlpha = .1;  % transparency of error bar shading

% reassign settings contained in opts
if exist('opts', 'var'); for i = 1:2:length(opts); s.(opts{i}) = opts{i+1}; end; end



% initializations
if ~iscell(eventTimes); eventTimes = {eventTimes}; end
plotEpochs = min(size(eventTimes{1}))==2;  % whehther this is a normal PSTH or interpolated epoch PSTH
numConditions = length(eventTimes);
if ischar(s.colors); s.colors = eval([s.colors '(numConditions)']); end % set colorspace if color is specified as a string

% load spike rates and times for all sessions (doing this first allows you to know the number of cells in advance, which helps with initializations)
temp = load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'neuralData.mat'), ...
    'timeStamps', 'spkRates');
timeStamps = temp.timeStamps;
spkRates = temp.spkRates(cellNum, :);
clear temp
fs = round(1/median(diff(timeStamps)));

% figure out x axis values
if ~plotEpochs
    times = s.xLims(1):(1/fs):s.xLims(2);
    s.xGridLength = length(times);
else
    % create one time axis for each condition
    times = cell(1,length(eventTimes));
    for i = 1:length(eventTimes)
        times{i} = linspace(0, median(diff(eventTimes{i},1,2)), s.xGridLength); % for this
    end
end

% get neural data for session
if s.normalize
    spkRates = zscore(spkRates);
    yLims = s.yLimsNormalized;
else
    rateMean = nanmean(spkRates);
    rateStd = nanstd(spkRates);
    yLims = [max(rateMean+s.yLimsNormalized(1)*rateStd,0), rateMean+s.yLimsNormalized(2)*rateStd];
end


cellData = cell(1, length(eventTimes));  % cell array that will contain one matrix of firing rates per condition
for eventType = 1:length(eventTimes)
    
    cellData{eventType} = nan(length(eventTimes{eventType}), s.xGridLength); % initialize data containers
    
    % get response for each trial of given eventType
    for event = 1:length(eventTimes{eventType})

        % get event responses
        if ~plotEpochs
            
            % get trial response
            startInd = find(timeStamps >= eventTimes{eventType}(event) + s.xLims(1), 1, 'first');
            if ~isempty(startInd)
                cellData{eventType}(event,:) = spkRates(startInd:startInd+length(times)-1);
            end

        % get epoch responses
        else
            % get trial response
            epochBins = timeStamps > eventTimes{eventType}(event,1) & ...
                        timeStamps < eventTimes{eventType}(event,2);
            if any(epochBins)
                x = timeStamps(epochBins);
                y = spkRates(epochBins);
                cellData{eventType}(event,:) = interp1(x, y, linspace(x(1), x(end), s.xGridLength));
            end
        end
    end
end



% plot!
hold on
if ~plotEpochs; x=times; else; x = mean(cat(1, times{:}),1); end % use average x axis across all conditions

for i = 1:numConditions 
    
    % plot individual trials
    if s.plotTrials>0
        randTrials = randperm(size(cellData{i},1), min(s.plotTrials, size(cellData{i},1)));
        for trial = randTrials
            plot(x, cellData{i}(trial,:), 'linewidth', 1, 'color', [s.colors(i,:) s.trialAlpha]);
        end
    end

    % plot mean and std
    if isequal(s.errorFcn, false) % plot only the mean
        plot(x, nanmean(cellData{i},1), 'linewidth', 3, 'color', s.colors(i,:));
    else
        shadedErrorBar(x, cellData{i}, {@nanmean, s.errorFcn}, ...
            'lineprops', {'linewidth', 3, 'color', s.colors(i,:)}, 'patchSaturation', s.errorAlpha);
    end
end

% pimp fig
set(gca, 'box', 'off', 'XLim', [x(1) x(end)]);
set(gca, 'YLim', yLims);

% add vertical line at time 0
if ~plotEpochs
    vertLine = line([0 0], get(gca, 'YLim'), 'color', [.6 .6 .6], 'lineWidth', 2);
    uistack(vertLine, 'bottom');
end

if numConditions>1 || s.showSampleSize
    sampleSizes = nan(1,numConditions);
    if isempty(s.conditionNames); s.conditionNames = cell(1,numConditions); end
    
    for i = 1:numConditions
        sampleSizes(i) = sum(~all(isnan(cellData{i}),2)); % number of trials where all of the entries in the firing rate matrix are not nan
        if s.showSampleSize; s.conditionNames{i} = sprintf('%s (%i)', s.conditionNames{i}, sampleSizes(i)); end
    end 
    
    for i = 1:numConditions; lines(i) = plot([nan nan], 'color', s.colors(i,:), 'LineWidth', 2); end % create dummy lines
    legend(lines, s.conditionNames, 'Box', 'off', 'Location', 'Best');
end
pause(.001)



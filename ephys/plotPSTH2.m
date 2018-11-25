function plotPSTH2(session, cellNum, eventTimes, conditionNames, conditionColors, plotStd)

% to do: interp all cells within a session at the same time (to speed
% things up) // plot multiple thing on top of eachother, which would let me
% get rid of bootstrappinfg too... // add vertical line option for epoch
% mode, eg to show swing and stance // accept previously computed data, and
% replot only, mwahahahaha // number plots // add buffer to begining and
% end of epoch plots

% !!! need to document


% settings
showFigure = false; % set to false if using this function to create subplots using separate code
yLimsNormalized = [-2 2];
yLims = [0 200];
stimPrePost = [-.1 .5];
fs = 1000;
normalize = false;
plotTrials = false;
trialsToPlot = 10;
plotControlDistribution = false;
controlColor = [0 0 0 .4];
xGridLength = 100; % number of points per trial for epoch interpolation

% initializations
if ~exist('conditionColors', 'var'); colors = hsv(length(eventTimes)); else; colors = conditionColors; end
if ~exist('plotStd', 'var'); plotStd = true; end
if ~iscell(eventTimes); eventTimes = {eventTimes}; end
plotEpochs = min(size(eventTimes{1}))==2;
numConditions = length(eventTimes);

if ~plotEpochs
    times = stimPrePost(1):(1/fs):stimPrePost(2);
    xGridLength = length(times);
else
    % create one time axis for each condition
    times = cell(1,length(eventTimes));
    for i = 1:length(eventTimes)
        times{i} = linspace(0, median(diff(eventTimes{i},1,2)), xGridLength);
    end
end


% load spike rates and times for all sessions (doing this first allows you to know the number of cells in advance, which helps with initializations)
temp = load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'neuralData.mat'), 'timeStamps', 'spkRates');
timeStamps = temp.timeStamps;
spkRates = temp.spkRates(cellNum, :);


% get neural data for session
if round(1/median(diff(timeStamps)))~=fs; disp('WARNING: expected frequency for instantaneous firing rate was not met!'); keyboard; end
if normalize; spkRates = zscore(spkRates); end
cellData = cell(1, length(eventTimes));
if plotControlDistribution; cellControlData = cell(1, length(eventTimes)); end

for eventType = 1:length(eventTimes)
    
    % initialize data containers
    eventData = nan(length(eventTimes{eventType}), xGridLength);
    if plotControlDistribution
        eventControlData = nan(length(eventTimes(eventType)), xGridLength);
    end
    
    % get response for each trial of given eventType
    for event = 1:length(eventTimes{eventType})

        % get event responses
        if ~plotEpochs
            % get trial response
            startInd = find(timeStamps >= eventTimes{eventType}(event) + stimPrePost(1), 1, 'first');
            if ~isempty(startInd)
                eventData(event,:) = spkRates(startInd:startInd+length(times)-1);
            end

            % get control distribution (bootstrap)
            if plotControlDistribution
                startInd = randi(length(timeStamps)-xGridLength);
                eventControlData(event,:) = spkRates(startInd:startInd+xGridLength-1);
            end

        % get epoch responses
        else
            % get trial response
            epochBins = timeStamps > eventTimes{eventType}(event,1) & ...
                        timeStamps < eventTimes{eventType}(event,2);
            if any(epochBins)
                x = timeStamps(epochBins);
                y = spkRates(epochBins);
                eventData(event,:) = interp1(x, y, linspace(x(1),x(end),xGridLength));
            end

            % get control distribution (bootstrap)
            if plotControlDistribution
                startTime = timeStamps(randi(length(timeStamps)-sum(epochBins)));
                epochBins = timeStamps > startTime & ...
                            timeStamps < startTime + diff(eventTimes{eventType}(event,:));
                if any(epochBins)
                    x = timeStamps(epochBins);
                    y = spkRates(epochBins);
                    eventControlData(event,:) = interp1(x, y, linspace(x(1),x(end),xGridLength));
                end
            end
        end
    end
    
    cellData{eventType} = eventData;
    if plotControlDistribution; cellControlData{eventType} = eventControlData; end
end



% plot!
% figure('Visible', showFigure, 'color', 'white', 'MenuBar', 'none', 'units', 'pixels', 'position', [2000 400 500 300]); hold on
% axes(parent)
if ~plotEpochs; x=times; else; x = mean(cat(1, times{:}),1); end % use average x axis across all conditions
meanTraces = cell(1, length(eventTimes));

% get control distribution
if plotControlDistribution
    
    controlData = cat(1, cellControlData{:});
    controlMean = nanmean(controlData, 1);
    controlStd = nanstd(controlData, 1);
    
    plot(x, controlMean, 'linewidth', 3, 'color', controlColor); hold on
    plot(x, controlMean + controlStd, 'linewidth', 0.5, 'color', controlColor)
    plot(x, controlMean - controlStd, 'linewidth', 0.5, 'color', controlColor)
end


for i = 1:numConditions
    
    condMean = nanmean(cellData{i},1);
    condStd = nanstd(cellData{i},1);
    
    % plot individual trials
    if plotTrials
        randTrials = randperm(size(cellData{i},1));
        randTrials = randTrials(1:min(trialsToPlot, size(cellData{i},1)));
        for trial = randTrials
            plot(x, cellData{i}(trial,:), 'linewidth', 1, 'color', [colors(i,:) 0.1]); hold on
        end
    end

    % plot mean and std
    meanTraces{i} = plot(x, condMean, 'linewidth', 3, 'color', colors(i,:)); hold on
    if plotStd
        plot(x, condMean + condStd, 'linewidth', 0.5, 'color', colors(i,:))
        plot(x, condMean - condStd, 'linewidth', 0.5, 'color', colors(i,:))
    end
end

% pimp fig
set(gca, 'box', 'off', 'XLim', [x(1) x(end)]);
if normalize; set(gca, 'YLim', yLimsNormalized); else; set(gca, 'YLim', yLims); end
if ~plotEpochs; line([0 0], get(gca, 'YLim'), 'color', controlColor); end
if numConditions>1; legend([meanTraces{:}], conditionNames, 'Box', 'off'); end



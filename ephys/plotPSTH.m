function allCellsData = plotPSTH(sessions, eventTimes, figTitle)

% to do: interp all cells within a session at the same time (to speed
% things up) // plot multiple thing on top of eachother, which would let me
% get rid of bootstrapping too... // add vertical line option for epoch
% mode, eg to show swing and stance // accept previously computed data, and
% replot only, mwahahahaha // number plots

% sessions: cell array of session names
% eventTimes: cell array with every entry a vector of event times in a
% session // if times are 2Xn matrices rather than vectors this function
% interpolates over epochs along a common x axis


% settings
yLimsNormalized = [-2 2];
yLims = [0 100];
stimPrePost = [-.5 2];
fs = 1000;
normalize = false;
sem = false;
plotTrials = false;
plotStd = true;
plotControlDistribution = true;
colors = hsv(length(sessions))*.8;
controlColor = [0 0 0 .4];
rows = 6;
xGridLength = 100; % number of points per trial for epoch interpolation


% initializations
plotEpochs = min(size(eventTimes{1}))==2;
if ~plotEpochs
    times = stimPrePost(1):(1/fs):stimPrePost(2);
    xGridLength = length(times);
else
    % create one time axis for 
    times = cell(1,length(sessions));
    for i = 1:length(sessions)
        times{i} = linspace(0, median(diff(eventTimes{i},1,2)), xGridLength);
    end
end
if plotControlDistribution; allSessionsControlData = cell(1,length(sessions)); end


% load spike rates and times for all sessions (doing this first allows you to know the number of cells in advance, which helps with initializations)
allSpkTimes = cell(1,length(sessions));
allSpkRates = cell(1,length(sessions));
for ses = 1:length(sessions)
    temp = load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{ses}, 'neuralData.mat'), 'spkTimes', 'spkRates');
    allSpkTimes{ses} = temp.spkTimes;
    allSpkRates{ses} = temp.spkRates;
end
totalCells = sum(cellfun(@(x) size(x,1), allSpkRates));
rows = min(totalCells, rows);
allCellsData = cell(1, totalCells);
plotControlDistribution; allCellsControlData = cell(1, totalCells);



% create PSTH matrices
disp('collecting data...')
cellInd = 1;
sessionIds = nan(1,totalCells);
for ses = 1:length(sessions)
    
    % get neural data for session
    if round(1/median(diff(allSpkTimes{ses}))) ~= fs; disp('WARNING: expected frequency for instantaneous firing rate was not met!'); keyboard; end
    if normalize; allSpkRates{ses} = zscore(allSpkRates{ses}); end
    
    % get events for each neuron
    for neuron = 1:size(allSpkRates{ses},1)
        fprintf('%.2f... ', cellInd/totalCells);
        cellData = nan(length(eventTimes{ses}), xGridLength); % event X time
        cellControlData = nan(length(eventTimes{ses}), xGridLength); % event X time
        
        % get spike rates for each event
        for event = 1:length(eventTimes{ses})
            
            % get event responses
            if ~plotEpochs
                % get trial response
                startInd = find(allSpkTimes{ses} >= eventTimes{ses}(event)+stimPrePost(1), 1, 'first');
                if ~isempty(startInd)
                    cellData(event,:) = allSpkRates{ses}(neuron,startInd:startInd+length(times)-1);
                end

                % get control distribution (bootstrap)
                if plotControlDistribution
                    startInd = randi(length(allSpkTimes{ses})-xGridLength);
                    cellControlData(event,:) = allSpkRates{ses}(neuron,startInd:startInd+xGridLength-1);
                end
                
            % get epoch responses
            else
                % get trial response
                epochBins = allSpkTimes{ses} > eventTimes{ses}(event,1) & ...
                            allSpkTimes{ses} < eventTimes{ses}(event,2);
                if any(epochBins)
                    x = allSpkTimes{ses}(epochBins);
                    y = allSpkRates{ses}(neuron, epochBins);
                    cellData(event,:) = interp1(x, y, linspace(x(1),x(end),xGridLength));
                end
                
                % get control distribution (bootstrap)
                if plotControlDistribution
                    startTime = allSpkTimes{ses}(randi(length(allSpkTimes{ses})-sum(epochBins)));
                    epochBins = allSpkTimes{ses} > startTime & ...
                                allSpkTimes{ses} < startTime + diff(eventTimes{ses}(event,:));
                    if any(epochBins)
                        x = allSpkTimes{ses}(epochBins);
                        y = allSpkRates{ses}(neuron, epochBins);
                        cellControlData(event,:) = interp1(x, y, linspace(x(1),x(end),xGridLength));
                    end
                end
            end
        end
        allCellsData{cellInd} = cellData;
        allCellsControlData{cellInd} = cellControlData;
        sessionIds(cellInd) = ses;
        cellInd = cellInd+1;
    end
end
fprintf('\n');


% compute means and standard deviations
cellAvgs = cellfun(@(x) squeeze(nanmean(x,1)), allCellsData, 'UniformOutput', 0); cellAvgs = cat(1, cellAvgs{:});
cellControlAvgs = cellfun(@(x) squeeze(nanmean(x,1)), allCellsControlData, 'UniformOutput', 0); cellControlAvgs = cat(1, cellControlAvgs{:});
cellStd = cellfun(@(x) squeeze(nanstd(x,1)), allCellsData, 'UniformOutput', 0); cellStd = cat(1, cellStd{:});
% cellStd = cellfun(@(x) squeeze(nanstd(x,1) / sqrt(size(x,1))), allCellsData, 'UniformOutput', 0); cellStd = cat(1, cellStd{:});
cellControlStd = cellfun(@(x) squeeze(nanstd(x,1)), allCellsControlData, 'UniformOutput', 0); cellControlStd = cat(1, cellControlStd{:});
% cellControlStd = cellfun(@(x) squeeze(nanstd(x,1) / sqrt(size(x,1))), allCellsControlData, 'UniformOutput', 0); cellControlStd = cat(1, cellControlStd{:});



% plot!
figure('name', figTitle, 'color', 'white', 'MenuBar', 'none', 'units', 'pixels', 'position', [2000 25 1000 800]);
cols = ceil(totalCells/rows);

for i = 1:totalCells
    
    subaxis(rows, cols, i, 'margin', .08, 'spacing', 0.05)
    [col,row] = ind2sub([cols,rows],i);
    lineColor = colors(sessionIds(i),:);
    if plotEpochs; x=times{sessionIds(i)}; else; x=times; end
    
    % plot individual trials
    if plotTrials
        for j = 1:size(allCellsData{i})
            plot(x, allCellsData{i}(j,:), 'linewidth', 1, 'color', [0 0 0 0.1]); hold on
        end
    end
    
    % plot bootstrap mean and std
    if plotControlDistribution
        plot(x, cellControlAvgs(i,:), 'linewidth', 3, 'color', controlColor); hold on
        plot(x, cellControlAvgs(i,:)+cellControlStd(i,:), 'linewidth', 0.5, 'color', controlColor)
        plot(x, cellControlAvgs(i,:)-cellControlStd(i,:), 'linewidth', 0.5, 'color', controlColor)
    end
    
    % plot mean and std
    plot(x, cellAvgs(i,:), 'linewidth', 3, 'color', lineColor); hold on
    if plotStd
        plot(x, cellAvgs(i,:)+cellStd(i,:), 'linewidth', 0.5, 'color', lineColor)
        plot(x, cellAvgs(i,:)-cellStd(i,:), 'linewidth', 0.5, 'color', lineColor)
    end
    
    
    set(gca, 'box', 'off', 'XLim', [x(1) x(end)]);
    if normalize; set(gca, 'YLim', yLimsNormalized); end % else; set(gca, 'YLim', yLims); end
    if (row~=rows || col~=1) && ~plotEpochs; set(gca, 'XColor', [1 1 1]); end
    line([0 0], get(gca, 'YLim'), 'color', controlColor)
end

disp('all done!')




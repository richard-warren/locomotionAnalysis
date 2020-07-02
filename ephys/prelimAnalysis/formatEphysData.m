function formatEphysData(session, varargin)
% work in progress // prepares spikes for matlab-land by getting spike
% times wrt spike clock, and maybe by computing instantaneous firing rates
% as well...


% settings
s.spkRateFs = 1000;     % sampling frequency of instantaneous firing rate
s.kernelRise = .005;  % rise and fall for double exponential kernel
s.kernelFall = .02;
s.kernelSig = .02;  % if a gaussian kernel is used
s.kernel = 'doubleExp';  % 'gauss', or 'doubleExp'
s.event = 'reward'; % what event is openEphys using for syncing
s.eventTimeName = 'rewardTimes'; % the name of paramters storing event(defined by s.event) time in spike clock


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs

addpath(fullfile(getenv('GITDIR'), 'analysis-tools'))
addpath(fullfile(getenv('GITDIR'), 'npy-matlab'))
files = dir(fullfile(getenv('OBSDATADIR'), 'sessions', session));
ephysFolder = files([files.isdir] & contains({files.name}, 'ephys_')).name;


% figure out which sync signal was used
ephysInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'ephysInfo.xlsx'), 'Sheet', 'ephysInfo');
s.event = ephysInfo.syncSignal{strcmp(session, ephysInfo.session)};
if strcmp(s.event, 'reward')
    s.eventTimeName = 'rewardTimes';
elseif strcmp(s.event, 'obsOn')
    s.eventTimeName = 'obsOnTimes';
else
    fprintf('WARNING: CANNOT figure out sync signal for openEphys and spike');
end
    

% get mapping from open ephys to spike times
[channel, eventTimes, info] = load_open_ephys_data_faster(...
    fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, 'all_channels.events'));
eventChannel = unique(channel); % assumes only one digital input is used!
eventEphysTimes = eventTimes(logical(info.eventId) & channel==eventChannel); % only take rising edge of event channel // !!! is the first variablee returned from load_open_ephys_data_faster really the identity of the event channel???
eventSpikeTimes = load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), s.eventTimeName); 
eventSpikeTimes = cell2mat(struct2cell(eventSpikeTimes));


if length(eventEphysTimes)~=length(eventSpikeTimes)
    fprintf('WARNING: %i events in spike and %i events in openEphys...', length(eventSpikeTimes), length(eventEphysTimes))


    if length(eventEphysTimes) > length(eventSpikeTimes)
        longString = eventEphysTimes;
        shortString = eventSpikeTimes;
    else
        longString = eventSpikeTimes;
        shortString = eventEphysTimes;
    end
    
    if longString(1) > shortString(1)
        sumDiff = [];
        matchPositionMatrix = nan(length(shortString), length(longString));
        shortStringMatrix = nan(length(shortString), length(longString));
        for i = 1:length(longString)
            difference = longString(i) - shortString(1);
            shortString = shortString + difference;
            shortStringMatrix(:, i) = shortString;
            matchPositionMatrix(:, i) = knnsearch(longString, shortString); % each colomn is one round of knnsearch
            sum = 0;
            for j = 1:length(shortString); sum = sum + abs(shortString(j) - longString(matchPositionMatrix(j, i))); end
            sumDiff(i) = sum;
        end
    else
        sumDiff = [];
        matchPositionMatrix = nan(length(shortString), length(longString));
        shortStringMatrix = nan(length(shortString), length(longString));
        for i = 1:length(longString)
            difference = shortString(1) - longString(i);
            shortString = shortString - difference;
            shortStringMatrix(:, i) = shortString;
            matchPositionMatrix(:, i) = knnsearch(longString, shortString); % each colomn is one round of knnsearch
            sum = 0;
            for j = 1:length(shortString); sum = sum + abs(shortString(j) - longString(matchPositionMatrix(j, i))); end
            sumDiff(i) = sum;
        end
    end
    
    % Get matched event time for both spike and ephys
    matchPosition = find(sumDiff == min(sumDiff));
    minTimeDiff = min(sumDiff);
    fprintf('\n Found match position in longer string, match position is %i ', matchPosition)
    fprintf('\n The minimal sum of time shifts for match position is %4.2f seconds\n', minTimeDiff)
    
    if length(eventEphysTimes) > length(eventSpikeTimes)
        matchEventSpikeTimes = eventSpikeTimes;
        matchEventEphysTimes = eventEphysTimes(matchPositionMatrix(:, matchPosition));
    else
        matchEventEphysTimes = eventEphysTimes;
        matchEventSpikeTimes = eventSpikeTimes(matchPositionMatrix(:, matchPosition));
    end
    
    
    % Validation Plot
    x1 = 1:length(longString);
    plot(x1, longString, '.', 'Color', [0.98 0.83 0.22], 'Markersize', 20)
    hold on
    x2 = 1:length(shortString);
    plot(matchPositionMatrix(:, matchPosition), shortStringMatrix(:, matchPosition), 'o', 'Color', [0.07 0.62 1], 'MarkerSize', 5)
    hold on
    plot(matchPositionMatrix(:, matchPosition), longString(matchPositionMatrix(:, matchPosition)), '.', 'Color', [1 0.36 0.40], 'MarkerSize', 10)
    legend('LongString-AllPoints', 'ShortString-AllPoints', 'LongString-MatchPoints');
    xlabel('Event Number');
    ylabel('Time');
    
else
    % [validOpenEBins, validSpikeBins] = deal(true(1,length(openEphysObsOnTimes)));
    disp('correct number of events detected!')
    matchEventSpikeTimes = eventSpikeTimes;
    matchEventEphysTimes = eventEphysTimes;
end


% Linear mapping from open ephys to spike event time
openEphysToSpikeMapping = polyfit(matchEventEphysTimes, matchEventSpikeTimes, 1); % get linear mapping from open ephys to spike

% check that predictions are accurate
predictedEventSpikeTimes = polyval(openEphysToSpikeMapping, matchEventEphysTimes);
if max(abs(predictedEventSpikeTimes - matchEventSpikeTimes)) > .001
    fprintf('OMG THERE WAS A PROBLEM MAPPING FROM OPEN EPHYS TO SPIKE TIMES!!! LOL\n');
else
    fprintf('MAPPING FROM OPEN EPHYS TO SPIKE TIMES SUCCESSFUL!!! LOL\n');
end
    


% get spike times for good units
[spkInds, unit_ids, ~] = getGoodSpkInds(session);
cellData = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'cellData.csv'));
if ~all(cellData.unit_id==unit_ids); disp('WARNING: callData.csv unit_ids do not match those in ephysFolder'); keyboard; end
ephysInfo = getSessionEphysInfo(session);

% restrict to good units
goodBins = logical([cellData.include]);
unit_ids = unit_ids(goodBins);
spkInds = spkInds(goodBins);
cellData = cellData(goodBins,:);

spkTimes = cell(1, length(unit_ids));
for i = 1:length(unit_ids)
    spkTimes{i}  = polyval(openEphysToSpikeMapping, ephysInfo.timeStamps(spkInds{i}));
end


% convert to instantaneous firing rate
minTime = min(cellfun(@(x) x(1), spkTimes)); % latest spike across all neurons
maxTime = max(cellfun(@(x) x(end), spkTimes)); % latest spike across all neurons

if strcmp(s.kernel, 'doubleExp')
    [~, timeStamps] = getFiringRateDoubleExp(spkTimes{1}, s.spkRateFs, s.kernelRise, s.kernelFall, [minTime maxTime]);
elseif strcmp(s.kernel, 'gauss')
    [~, timeStamps] = getFiringRateGaussian(spkTimes{1}, s.spkRateFs, s.kernelSig, [minTime maxTime]);
end
spkRates = nan(length(spkTimes), length(timeStamps));

for i = 1:length(spkTimes)
    
    if strcmp(s.kernel, 'doubleExp')
        spkRates(i,:) = getFiringRateDoubleExp(spkTimes{i}, s.spkRateFs, s.kernelRise, s.kernelFall, [minTime maxTime]);
    elseif strcmp(kernel, 'gauss')
        spkRates(i,:) = getFiringRateGaussian(spkTimes{i}, s.spkRateFs, s.kernelSig, [minTime maxTime]);
    end
    
    % get min and max time for cell
    cellMinTime = polyval(openEphysToSpikeMapping, cellData.timeStart(i)*60);
    if strcmp(cellData.timeEnd(i), 'max')
        cellMaxTime = timeStamps(end);
    else
        if iscell(cellData.timeEnd(i)); timeEnd = str2double(cellData.timeEnd(i)); else; timeEnd = cellData.timeEnd(i); end
        cellMaxTime = polyval(openEphysToSpikeMapping, timeEnd*60);
    end
    
    % remove spikes that are out of min and max times
    spkRates(i, timeStamps<cellMinTime | timeStamps>cellMaxTime) = nan;
    spkTimes{i} = spkTimes{i}(spkTimes{i}>cellMinTime & spkTimes{i}<cellMaxTime);
    
end


save(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'neuralData.mat'), ...
     'spkRates', 'spkTimes', 'timeStamps', 'unit_ids', 'openEphysToSpikeMapping')
disp('all done!')





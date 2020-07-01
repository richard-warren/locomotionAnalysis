function formatEphysData(session, varargin)
% work in progress // prepares spikes for matlab-land by getting spike
% times wrt spike clock, and maybe by computing instantaneous firing rates
% as well...


% settings
s.spkRateFs = 1000;      % sampling frequency of instantaneous firing rate
s.kernelRise = .005;     % (s) rise for double exponential kernel
s.kernelFall = .02;      % (s) fall for double exponential kernel
s.kernelSig = .02;       % (s) if a gaussian kernel is used
s.kernel = 'doubleExp';  % 'gauss', or 'doubleExp'


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs
addpath(fullfile(getenv('GITDIR'), 'analysis-tools'))
addpath(fullfile(getenv('GITDIR'), 'npy-matlab'))
ephysFolder = dir(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'ephys_*'));
ephysFolder = ephysFolder.name;
ephysInfo = getSessionEphysInfo(session);

% get mapping from open ephys to spike times
[channel, eventTimes, info] = load_open_ephys_data_faster(...
    fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, 'all_channels.events'));  % load event data
eventChannel = unique(channel);  % assumes only one digital input is used!
syncEphysTimes = eventTimes(logical(info.eventId) & channel==eventChannel); % only take rising edge of event channel
syncSpikeTimes = load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'run.mat'), ephysInfo.syncSignal);
syncSpikeTimes = syncSpikeTimes.(ephysInfo.syncSignal).times(syncSpikeTimes.(ephysInfo.syncSignal).level==1);


if length(syncEphysTimes)~=length(syncSpikeTimes)
    fprintf('%s: WARNING! %i events in spike and %i events in openEphys...\n', ...
        session, length(syncSpikeTimes), length(syncEphysTimes))
    
    % find signal offset
    % ------------------
    
    % make kernel
    fs = 1000;
    sig = .1;  % (seconds) std for gaussian
    kernel = arrayfun(@(x) (1/(sig*sqrt(2*pi))) * exp(-.5*(x/sig)^2), (-sig*10):(1/fs):(sig*10));
    
    % convolve
    tLims = [syncSpikeTimes(1)-range(syncEphysTimes) syncSpikeTimes(end)+range(syncEphysTimes)];
    t = tLims(1) : 1/fs : tLims(2);
    spikeBins = histcounts(syncSpikeTimes, t);
    spikeConv = conv(spikeBins, kernel, 'same');
    ephysBins = histcounts(syncEphysTimes, t);
    ephysConv = conv(ephysBins, kernel, 'same');
    
    % align signals
    [r, lags] = xcorr(spikeConv, ephysConv);
    [~, maxInd] = max(r);
    lag = lags(maxInd)/fs;
    syncEphysTemp = syncEphysTimes + lag;
    syncSpikeTemp = syncSpikeTimes;
    
    
    % find matching events
    % --------------------
    maxDiff = .1;
    matchedInds = nan(2,0);  % first row spike, second row ephys
    
    while ~all(isnan(syncEphysTemp)) && ~all(isnan(syncSpikeTemp))
        
        diffs = abs(syncSpikeTimes - syncEphysTemp');  % abs(diffs) between each pair of events (spkInds X ephysInds)
        [minDiff, minInd] = min(diffs, [], 'all', 'linear');
        if minDiff<maxDiff
            [spikeInd, ephysInd] = ind2sub(size(diffs), minInd);
            syncSpikeTemp(spikeInd) = nan;
            syncEphysTemp(ephysInd) = nan;
            matchedInds(:,end+1) = [spikeInd; ephysInd];
        else
            break
        end
    end
    
    
    % plot
    % ----
    figure('name', sprintf('%s: spike / openephys event alignment', session), ...
        'color', 'white', 'position', [90.00 654.00 1769.00 225.00]); hold on
    plot([syncSpikeTimes(matchedInds(1,:)), syncEphysTimes(matchedInds(2,:))+lag], ...
        [1 0], 'LineWidth', 1, 'color', [0 0 0 .4])  % lines connecting matched events
    scatter(syncSpikeTimes(matchedInds(1,:)), ones(1,size(matchedInds,2)), 50, [0 0.44 0.74], 'filled')
    scatter(syncEphysTimes(matchedInds(2,:))+lag, zeros(1,size(matchedInds,2)), 50, [.85 .32 .10], 'filled')
    
    spkUnmatchedInds = find(~ismember(1:length(syncSpikeTimes), matchedInds(1,:)));
    scatter(syncSpikeTimes(spkUnmatchedInds), ones(1,length(spkUnmatchedInds)), 100, 'red')
    ephysUnmatchedInds = find(~ismember(1:length(syncEphysTimes), matchedInds(2,:)));
    scatter(syncEphysTimes(ephysUnmatchedInds)+lag, zeros(1,length(ephysUnmatchedInds)), 100, 'red')
    
    set(gca, 'ytick', [0 1], 'YTickLabel', {'open ephys', 'spike'}, 'ylim', [-1 2])
    xlabel('spike times (s)')
    pause(.1)
    
    
    % update times
    syncSpikeTimes = syncSpikeTimes(matchedInds(1,:));
    syncEphysTimes = syncEphysTimes(matchedInds(2,:));
end


% linear mapping from open ephys to spike event time
% todo: the mapping should maybe be piecewise, eg using interp1()
openEphysToSpikeMapping = polyfit(syncEphysTimes, syncSpikeTimes, 1); % get linear mapping from open ephys to spike

% check that predictions are accurate
predictedEventSpikeTimes = polyval(openEphysToSpikeMapping, syncEphysTimes);
if max(abs(predictedEventSpikeTimes - syncSpikeTimes)) > .001
    fprintf('OMG THERE WAS A PROBLEM MAPPING FROM OPEN EPHYS TO SPIKE TIMES!!! LOL\n');
end
    

% get spike times for good units
[spkInds, unit_ids] = getGoodSpkInds(session);
cellData = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'cellData.csv'));
if ~all(cellData.unit_id==unit_ids); disp('WARNING: callData.csv unit_ids do not match those in ephysFolder'); keyboard; end


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

settings = s;
save(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'neuralData.mat'), ...
     'spkRates', 'spkTimes', 'timeStamps', 'unit_ids', 'openEphysToSpikeMapping', 'settings')
disp('all done!')





function formatEphysData(session)
% work in progress // prepares spikes for matlab-land by getting spike
% times wrt spike clock, and maybe by computing instantaneous firing rates
% as well...


% settings
spkRateFs = 1000;     % sampling frequency of instantaneous firing rate
spkRateKernSig = .02; % kernel used to determine instantaneous firing rate
% obsOnChannel = 0;

% initializations
addpath(fullfile(getenv('GITDIR'), 'analysis-tools'))
addpath(fullfile(getenv('GITDIR'), 'npy-matlab'))
files = dir(fullfile(getenv('OBSDATADIR'), 'sessions', session));
ephysFolder = files([files.isdir] & contains({files.name}, 'ephys_')).name;


% get mapping from open ephys to spike times
[channel, openEphysObsOnTimes, info] = load_open_ephys_data_faster(...
    fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, 'all_channels.events'));
obsOnChannel = unique(channel); % assumes only one digital input is used!
openEphysObsOnTimes = openEphysObsOnTimes(logical(info.eventId) & channel==obsOnChannel); % only take rising edge of event channel // !!! is the first variablee returned from load_open_ephys_data_faster really the identity of the event channel???
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), 'obsOnTimes');


if length(openEphysObsOnTimes)~=length(obsOnTimes)
    fprintf('WARNING: %i obsOnTimes in spike and %i obsOnTimes in openEphys...', length(obsOnTimes), length(openEphysObsOnTimes))
    
    validOpenEBins = false(1,length(openEphysObsOnTimes));
    validSpikeBins = false(1, length(obsOnTimes));
    openEphysObsOnTimesShifted = openEphysObsOnTimes-openEphysObsOnTimes(1)+obsOnTimes(1); % assumes the first event is shared for open ephys and spike
    for i = 1:length(obsOnTimes)
        [dt, closestBin] = min(abs(openEphysObsOnTimesShifted - obsOnTimes(i)));
        if dt<.1
            validOpenEBins(closestBin) = true;
            validSpikeBins(i) = true;
        end
    end
    
    fprintf(' but fixed using a hacky hack\n')
else
    [validOpenEBins, validSpikeBins] = deal(true(1,length(openEphysObsOnTimes)));
    disp('correct number of events detected!')
end
% fprintf('\n')


openEphysToSpikeMapping = polyfit(openEphysObsOnTimes(validOpenEBins), obsOnTimes(validSpikeBins), 1); % get linear mapping from open ephys to spike

% check that predictions are accurate
predictedObsOnTimes = polyval(openEphysToSpikeMapping, openEphysObsOnTimes(validOpenEBins));
if max(abs(predictedObsOnTimes - obsOnTimes)) > .001
    fprintf('OMG THERE WAS A PROBLEM MAPPING FROM OPEN EPHYS TO SPIKE TIMES!!! LOL\n');
    keyboard
end




% get spike times for good units
[spkInds, unit_ids] = getGoodSpkInds(session);
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

[~, timeStamps] = getFiringRate(spkTimes{1}, spkRateFs, spkRateKernSig, [minTime maxTime]);
spkRates = nan(length(spkTimes), length(timeStamps));

for i = 1:length(spkTimes)
    spkRates(i,:) = getFiringRate(spkTimes{i}, spkRateFs, spkRateKernSig, [minTime maxTime]);
    
    % get min and max time for cell
    cellMinTime = polyval(openEphysToSpikeMapping, cellData.timeStart(i)*60);
    if strcmp(cellData.timeEnd(i), 'max')
        cellMaxTime = timeStamps(end);
    else
        cellMaxTime = polyval(openEphysToSpikeMapping, str2num(cellData.timeEnd{i})*60);
    end
    
    % remove spikes that are out of min and max times
    spkRates(i, timeStamps<cellMinTime | timeStamps>cellMaxTime) = nan;
    spkTimes{i} = spkTimes{i}(spkTimes{i}>cellMinTime & spkTimes{i}<cellMaxTime);
    
end

save(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'neuralData.mat'), ...
    'spkRates', 'spkTimes', 'timeStamps', 'unit_ids')
disp('all done!')




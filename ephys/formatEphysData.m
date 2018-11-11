function formatEphysData(session)
% work in progress // prepares spikes for matlab-land by getting spike
% times wrt spike clock, and maybe by computing instantaneous firing rates
% as well...


% settings
spkRateFs = 1000;     % sampling frequency of instantaneous firing rate
spkRateKernSig = .02; % kernel used to determine instantaneous firing rate
obsOnChannel = 0;

% initializations
addpath(fullfile(getenv('GITDIR'), 'analysis-tools'))
addpath(fullfile(getenv('GITDIR'), 'npy-matlab'))
files = dir(fullfile(getenv('OBSDATADIR'), 'sessions', session));
ephysFolder = files([files.isdir] & contains({files.name}, 'ephys_')).name;


% get mapping from open ephys to spike times
[channel, openEphysObsOnTimes, info] = load_open_ephys_data_faster(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, 'all_channels.events'));
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
    
    fprintf(' but fixed using a hacky hack')
else
    [validOpenEBins, validSpikeBins] = deal(true(1,length(openEphysObsOnTimes)));
    disp('correct number of events detected!')
end
fprintf('\n')


openEphysToSpikeMapping = polyfit(openEphysObsOnTimes(validOpenEBins), obsOnTimes(validSpikeBins), 1); % get linear mapping from open ephys to spike

% check that predictions are accurate
predictedObsOnTimes = polyval(openEphysToSpikeMapping, openEphysObsOnTimes(validOpenEBins));
if max(abs(predictedObsOnTimes - obsOnTimes)) > .001
    fprintf('OMG THERE WAS A PROBLEM MAPPING FROM OPEN EPHYS TO SPIKE TIMES!!! LOL\n');
    keyboard
end




% get spike times for good units
spkTimes = readNPY(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, 'spike_times.npy'));
clusters = readNPY(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, 'spike_clusters.npy'));
clusterGroups = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, 'cluster_groups.csv'));
unit_ids = clusterGroups.cluster_id(strcmp(clusterGroups.group, 'good'));

unitTimes = cell(1, length(unit_ids));
[~, timeStamps] = load_open_ephys_data_faster(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, '107_CH1.continuous')); % assumes timestamps for all channels are the same... this is an okay assumption, right???
for i = 1:length(unit_ids)
    unitSmps = spkTimes(clusters==unit_ids(i));
    unitTimes{i}  = polyval(openEphysToSpikeMapping, timeStamps(unitSmps));
end


% convert to instantaneous firing rate

minTime = min(cellfun(@(x) x(1), unitTimes)); % latest spike across all neurons
maxTime = max(cellfun(@(x) x(end), unitTimes)); % latest spike across all neurons

[spkRatesTemp, spkTimes] = getFiringRate(unitTimes{1}, spkRateFs, spkRateKernSig, [minTime maxTime]);
spkRates = nan(length(unitTimes), length(spkTimes));
spkRates(1,:) = spkRatesTemp;
for i = 2:length(unitTimes)
    spkRates(i,:) = getFiringRate(unitTimes{i}, spkRateFs, spkRateKernSig, [minTime maxTime]);
end

save(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'neuralData.mat'), ...
    'spkRates', 'spkTimes', 'unitTimes')
% rmpath(fullfile(getenv('GITDIR'), 'analysis-tools'))
% rmpath(fullfile(getenv('GITDIR'), 'npy-matlab'))












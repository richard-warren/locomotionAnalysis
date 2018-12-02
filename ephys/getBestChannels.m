function [bestChannels, unit_ids] = getBestChannels(session, ephysInfo)

% takes a sample of all spikes, computes median waveform, and picks
% channels with best peak to peak amplitude

% settings
spkWindow = [-.25 .5]; % ms pre and post spike time to plot
spkNum = 500; % take this many of all spikes to analyze (to save system resources)



% initializations
if ~exist('ephysInfo', 'var'); ephysInfo = getSessionEphysInfo(session); end
spkWindowInds = int64((spkWindow(1)/1000*ephysInfo.fs) : (spkWindow(2)/1000*ephysInfo.fs));

% function to extract voltage from binary file, and load data
getVoltage = @(data, channel, inds) double(data.Data.Data(channel,inds))*ephysInfo.bitVolts;
data = memmapfile(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysInfo.ephysFolder, [ephysInfo.fileNameBase '_CHs.dat']), ...
    'Format', {'int16', [ephysInfo.channelNum, ephysInfo.smps], 'Data'}, 'Writable', false);

[allSpkInds, unit_ids] = getGoodSpkInds(session); % get spike times for good units

% get channel map
load(fullfile(getenv('OBSDATADIR'), 'ephys', 'channelMaps', 'kiloSort', [ephysInfo.mapFile '.mat']), 'connected');
connectedChannels = find(connected);
bestChannels = nan(1,length(unit_ids));

for c = 1:length(unit_ids)
    
    % extract waveform across all connected channels
    spkIndsSub = allSpkInds{c}(round(linspace(2, length(allSpkInds{c}), spkNum))); % get subpopulation of spikes evenly spaced out
    spkIndsSubAll = uint64(repmat(spkWindowInds,1,length(spkIndsSub)) + int64(repelem(spkIndsSub, length(spkWindowInds)))'); % inds for all all smps within all spikes
    allWaveforms = getVoltage(data, connectedChannels, spkIndsSubAll);
    
    % reshape
    allWaveforms = reshape(allWaveforms, length(connectedChannels), length(spkWindowInds), []);
    allWaveforms = permute(allWaveforms, [3 1 2]);
    
    % get best channel
    medianWaveform = squeeze(mean(allWaveforms,1));
    [~, maxInd] = max(peak2peak(medianWaveform,2));
    bestChannels(c) = connectedChannels(maxInd);
end





function [bestChannels, unit_ids] = getBestChannels(session)

% gets best channels for each unit by loading it from 'cluster_info.tsv' if
% that file exists (new sessions), or by computing from scratch (old
% sessions) by taking a sample of all spikes, computing median waveform,
% and picking channel with highest peak to peak amplitude

% settings
spkWindow = [-.25 .5]; % ms pre and post spike time to plot
spkSamples = 500; % take this many of all spikes to analyze (to save system resources)


% initializations
ephysInfo = getSessionEphysInfo(session);
ephysFolder = fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysInfo.ephysFolder);
spkWindowInds = int64((spkWindow(1)/1000*ephysInfo.fs) : (spkWindow(2)/1000*ephysInfo.fs));
loadBestChannels = exist(fullfile(ephysFolder, 'cluster_group.tsv'), 'file');


if loadBestChannels
    clusterInfo = tdfread(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysInfo.ephysFolder, 'cluster_info.tsv'));
    bins = all(clusterInfo.group(:,1:4)=='good',2);  % this is a bit of a hack // should really convert all rows to strings and use strcmp()
    unit_ids = clusterInfo.id(bins);
    
    if ismember('channel', fieldnames(clusterInfo))  % kilosort2
        bestChannels = clusterInfo.channel(bins) + 1;
    elseif ismember('ch', fieldnames(clusterInfo))  % kilosort1
        bestChannels = clusterInfo.ch(bins) + 1;
    end
else

    % function to extract voltage from binary file, and load data
    getVoltage = @(data, channel, inds) double(data.Data.Data(channel,inds))*ephysInfo.bitVolts;
    data = memmapfile(fullfile(ephysFolder, [ephysInfo.fileNameBase '_CHs.dat']), ...
        'Format', {'int16', [ephysInfo.channelNum, ephysInfo.smps], 'Data'}, 'Writable', false);

    [allSpkInds, unit_ids] = getGoodSpkInds(session); % get spike times for good units

    % get channel map
    load(fullfile(getenv('OBSDATADIR'), 'ephys', 'channelMaps', 'kiloSort', [ephysInfo.mapFile '.mat']), 'connected');
    connectedChannels = find(connected);
    bestChannels = nan(1,length(unit_ids));

    for c = 1:length(unit_ids)

        % extract waveform across all connected channels
        spkIndsSub = allSpkInds{c}(round(linspace(2, length(allSpkInds{c}), spkSamples))); % get subpopulation of spikes evenly spaced out
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
end




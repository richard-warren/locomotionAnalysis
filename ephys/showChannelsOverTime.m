function showChannelsOverTime(session, timeBinNum, showSortedSpikes, figureName, bestChannels)

% given session name, shows traces of all channels on a probe as rows in a
% large column // does this for evenly spaced intervals over the duration
% of recording // use to assess drift!

% TO DO: show only best channel spikes


% settings
% highPassFreq = 300;
yLims = [-400 400]; % microvols
if ~exist('timeBinNum', 'var'); timeBinNum = 4; end
if ~exist('showSortedSpikes', 'var'); showSortedSpikes = false; end
windowSize = .2/timeBinNum; % length of window in which to show spikes
if ~exist('figureName', 'var'); figureName = fullfile(getenv('OBSDATADIR'), 'figures', 'ephys', 'drift', [session '.png']); end
spkWindow = [-.5 1]; % ms pre and post spike time to plot


% initializations
ephysInfo = getSessionEphysInfo(session);
for i = fieldnames(ephysInfo)'; eval([i{1} '=ephysInfo.' i{1} ';']); end % extract field names as variables
load(fullfile(getenv('OBSDATADIR'), 'ephys', 'channelMaps', 'kilosort', [mapFile '.mat']), ...
    'xcoords', 'ycoords', 'connected')
[~, sortIndsTemp] = sort(ycoords); sortInds = nan(1,channelNum);
for i = 1:channelNum; sortInds(i) = find(sortIndsTemp==i); end % a hack, this works but dont really understand why
    
% timesSub = linspace(0,windowSize/fs,windowSize*fs);
if showSortedSpikes
    [spkInds, unit_ids] = getGoodSpkInds(session);
    spkWindowInds = int64((spkWindow(1)/1000*fs) : (spkWindow(2)/1000*fs));
    colors = hsv(length(spkInds));
else
    colors = hsv(channelNum);
end

% function to extract voltage from binary file
getVoltage = @(data, channel, inds) ...
    double(data.Data.Data(channel,inds))*bitVolts; % extract voltage from memmapfile, converting to votlage, highpassing, and only return specific channel

% load data
data = memmapfile(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, [fileNameBase '_CHs.dat']), ...
    'Format', {'int16', [channelNum, smps], 'Data'}, 'Writable', false);


% plot all channels at dft time intervals
figure('color', 'white', 'Units', 'pixels', 'Position', [2000 20 1800 950]);
timeBinEdges = linspace(min(timeStamps), max(timeStamps), timeBinNum+1);
timeStarts = timeBinEdges(1:end-1);
traces = nan(channelNum, length(timeStarts), windowSize*fs);

% collect traces)
for j = 1:length(timeStarts)
    startInd = find(timeStamps>timeStarts(j),1,'first');
    inds = startInd:startInd+windowSize*fs-1;
    traces(:,j,:) = getVoltage(data, 1:channelNum, inds);
end


% plot traces
timesSub = linspace(0, windowSize, windowSize*fs);
offsets = sortInds*(range(yLims));
for i = 1:channelNum
    for j = 1:length(timeStarts)
        subplot(1,length(timeStarts),j); hold on
        if showSortedSpikes; color=[.5 .5 .5]; else; color=colors(i,:); end
        
        plot(timesSub, squeeze(traces(i,j,:)) + offsets(i), 'Color', color);
        
        if connected(i); color='black'; else; color='red'; end
        text(timesSub(1), offsets(i), num2str(i), 'Color', color)
    end
end


% plot spikes on top of traces (if showSortedSpikes)
if showSortedSpikes
    lines = nan(1, length(unit_ids));
    for i = 1:length(unit_ids)
        for j = 1:length(timeStarts)
            subplot(1,length(timeStarts),j); hold on
            startInd = int64(find(timeStamps>timeStarts(j),1,'first'));

            % overlay spikes
            traceSpkInds = spkInds{i}(timeStamps(spkInds{i})>=timeBinEdges(j) & timeStamps(spkInds{i})<timeBinEdges(j)+windowSize);
            spkBins = repmat(spkWindowInds,length(traceSpkInds),1) + int64(traceSpkInds); % matrix where each row is inds for given spike
            spkBins = reshape(spkBins',[],1);
            spkBins = spkBins-startInd;
            spkBins = spkBins(spkBins<size(traces,3) & spkBins>0);

            if exist('bestChannels', 'var'); channelsToShow=bestChannels(i); else; channelsToShow=1:channelNum; end
            for k = channelsToShow
                spikesTrace = nan(1,size(traces,3));

                spikesTrace(spkBins) = traces(k,j,spkBins);
                lines(i) = plot(timesSub, spikesTrace + offsets(k), ...
                    'Color', [colors(i,:) .6], 'LineWidth', 1.5); % overlay spikes! lol
            end
        end
    end
end


% pimp fig
for i = 1:length(timeStarts)
    subplot(1,length(timeStarts),i);
    
    set(gca, 'XLim', [0 timesSub(end)], ...
        'YLim', [0 (channelNum+.5)*range(yLims)], ...
        'Visible', 'off')
    
    text(range(timesSub)*.1, (channelNum+1)*(range(yLims)), ...
        [num2str(round((timeStarts(i)-min(timeStarts))/60)) ' minutes'], ...
        'FontWeight', 'bold')
end
if showSortedSpikes; legend(lines, cellstr(num2str(unit_ids)), 'location', 'northeast'); end


% save that ish
saveas(gcf, figureName);

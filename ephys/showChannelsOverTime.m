function showChannelsOverTime(session, timeBinNum, showSortedSpikes, figureName, bestChannels)

% given session name, shows traces of all channels on a probe as rows in a
% large column // does this for evenly spaced intervals over the duration
% of recording // use to assess drift!

% TO DO: show only best channel spikes


% settings
% highPassFreq = 300;
yLims = [-500 500]; % microvols
if ~exist('timeBinNum', 'var'); timeBinNum = 4; end
if ~exist('showSortedSpikes', 'var'); showSortedSpikes = false; end
windowSize = 2.0/timeBinNum; % the length of the time window for plotting traces. eg, for each time interval, plot 0.5 sec of the data
if ~exist('figureName', 'var'); figureName = fullfile(getenv('OBSDATADIR'), 'figures', 'ephys', 'drift', [session '.png']); end
spkWindow = [-.5 1]; % ms pre and post spike time to plot


% initializations
ephysInfo = getSessionEphysInfo(session);
for i = fieldnames(ephysInfo)'; eval([i{1} '=ephysInfo.' i{1} ';']); end % extract field names as variables
load(fullfile(getenv('OBSDATADIR'), 'ephys', 'channelMaps', 'kilosort', [mapFile '.mat']), ...
    'xcoords', 'ycoords', 'connected', 'channelNum_OpenEphys')



% transfer the openephys channel order into its physical location order on the probe
% in order to determin the spacing offset for every channel for the drift plots.
sortedInds = [];
for i = 1:64
    temp = find(channelNum_OpenEphys == i); 
    sortedInds = [sortedInds; 65 - temp];
end
    


if showSortedSpikes
    [spkInds, unit_ids] = getGoodSpkInds(session);
    spkWindowInds = int64((spkWindow(1)/1000*fs) : (spkWindow(2)/1000*fs));
    colors = hsv(length(spkInds));
else
    colors = parula(channelNum);
end



% function to extract voltage from binary file
getVoltage = @(data, channel, inds) ...
    double(data.Data.Data(channel,inds))*bitVolts; % extract voltage from memmapfile, converting to votlage, highpassing, and only return specific channel

% load data
data = memmapfile(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, [fileNameBase '_CHs.dat']), ...
    'Format', {'int16', [channelNum, smps], 'Data'}, 'Writable', false);


% plot all channels at dft time intervals
figure('color', 'white', 'Units', 'pixels', 'position', get(0,'ScreenSize'));

timeBinEdges = linspace(min(timeStamps), max(timeStamps), timeBinNum+1);
timeStarts = timeBinEdges(1:end-1);
traces = nan(channelNum, length(timeStarts), windowSize*fs);


% collect traces
for j = 1:length(timeStarts)
    startInd = find(timeStamps>timeStarts(j),1,'first');
    inds = startInd:startInd+windowSize*fs-1;
    traces(:,j,:) = getVoltage(data, 1:channelNum, inds);
end


% plot traces
timesSub = linspace(0, windowSize, windowSize*fs); % xlim for plotting the traces
offsets = sortedInds*(range(yLims)); % set offsets for every channel based on its physical location on the probe
for i = 1:channelNum
    for j = 1:length(timeStarts)
        subplot(1,length(timeStarts),j); hold on
        
        hold on;
        if showSortedSpikes; color=[.5 .5 .5]; else; color=colors(i,:); end
        
        plot(timesSub, squeeze(traces(i,j,:)) + offsets(i), 'Color', color);
        
        if connected(i); color='black'; else; color='red'; end
        hold on
        text(timesSub(1), offsets(i), ['(' num2str(i-1) ')' num2str(65 - offsets(i)/1000)], 'Color', color)
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

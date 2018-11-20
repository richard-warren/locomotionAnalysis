function showChannelsOverTime(session)

% given session name, shows traces of all channels on a probe as rows in a
% large column // does this for evenly spaced intervals over the duration
% of recording // use to assess drift!


% settings
% highPassFreq = 300;
yLims = [-400 400]; % microvols
windowSize = .5; % length of window in which to show spikes
timeBins = 5;


% initializations
addpath(fullfile(getenv('GITDIR'), 'analysis-tools'))
files = dir(fullfile(getenv('OBSDATADIR'), 'sessions', session));
ephysFolder = files([files.isdir] & contains({files.name}, 'ephys_')).name;
contFiles = dir(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, '*.continuous'));
channelNum = length(contFiles);
fileNameBase = contFiles(1).name(1:3);

% get channel mapping
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
ephysInfo = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', 'ephysInfo.xlsx'), 'Sheet', 'ephysInfo');
warning('on', 'MATLAB:table:ModifiedAndSavedVarnames')
mapFile = ephysInfo.map{strcmp(session, ephysInfo.session)};
load(fullfile(getenv('OBSDATADIR'), 'ephys', 'channelMaps', 'kilosort', [mapFile '.mat']), ...
    'xcoords', 'ycoords', 'connected')
[~, sortInds] = sort(ycoords);


% get fs, microvolts conversion factor, and number of samples
[~, times, info] = load_open_ephys_data_faster(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, [fileNameBase '_CH1.continuous']));
fs = info.header.sampleRate;
bitVolts = info.header.bitVolts;
file = fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, [fileNameBase '_CHs.dat']);
temp = dir(file);
smps = temp.bytes/2/channelNum; % 2 bytes per sample
timesSub = linspace(0,windowSize/fs,windowSize*fs);

% function to extract voltage from binary file
getVoltage = @(data, channel, inds) ...
    double(data.Data.Data(channel,inds))*bitVolts; % extract voltage from memmapfile, converting to votlage, highpassing, and only return specific channel

% load data
data = memmapfile(file, 'Format', {'int16', [channelNum, smps], 'Data'}, 'Writable', false);


% plot all channels at dft time intervals
figure('color', 'white', 'Units', 'pixels', 'Position', [2000 20 1800 950]);
timeBinEdges = linspace(min(times), max(times), timeBins+1);
timeStarts = timeBinEdges(1:end-1);
traces = nan(channelNum, length(timeStarts), windowSize*fs);

% collect traces)
for j = 1:length(timeStarts)
    startInd = find(times>timeStarts(j));
    inds = startInd:startInd+windowSize*fs-1;
    traces(:,j,:) = getVoltage(data, 1:channelNum, inds);
end


% plot traces
for i = 1:channelNum
    for j = 1:length(timeStarts)
        subplot(1,length(timeStarts),j); hold on
        plot(timesSub, squeeze(traces(i,j,:)) + find(i==sortInds)*(range(yLims)));
        if connected(i); color='black'; else; color='red'; end
        text(timesSub(1), find(i==sortInds)*(range(yLims)), num2str(i), 'Color', color)
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

% save that ish
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures', 'ephys', 'drift', [session '.png']));

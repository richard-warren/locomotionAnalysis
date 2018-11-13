function showChannelsOverTime(session, minuteIntervals)

% given session name, shows traces of all channels on a probe as rows in a
% large column // does this for evenly spaced intervals over the duration
% of recording // use to assess drift!

% !!! fix dorsovental sequencing!


% settings
highPassFreq = 300;
yLims = [-400 400]; % microvols
windowSize = .5; % length of window in which to show spikes


% initializations
timeIntervals = minuteIntervals*60; % convert to seconds
disp('initializing data...')
files = dir(fullfile(getenv('OBSDATADIR'), 'sessions', session));
ephysFolder = files([files.isdir] & contains({files.name}, 'ephys_')).name;
contFiles = dir(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, '*.continuous'));
channelNum = length(contFiles);
fileNameBase = contFiles(1).name(1:3);

% get channel mapping
ephysInfo = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', 'ephysInfo.xlsx'), 'Sheet', 'ephysInfo');
mapFile = ephysInfo.map{strcmp(session, ephysInfo.session)};
load(fullfile(getenv('OBSDATADIR'), 'ephys', 'channelMaps', 'kilosort', [mapFile '.mat']), ...
    'xcoords', 'ycoords')
[~, sortInds] = sort(ycoords);


% get fs, microvolts conversion factor, and number of samples
[~, ~, info] = load_open_ephys_data_faster(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, [fileNameBase '_CH1.continuous']));
fs = info.header.sampleRate;
bitVolts = info.header.bitVolts;
file = fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, [fileNameBase '_CHs.dat']);
temp = dir(file);
smps = temp.bytes/2/channelNum; % 2 bytes per sample
times = linspace(0,smps/fs,smps);
timesSub = linspace(0,windowSize/fs,windowSize*fs);

% function to extract voltage from binary file
getVoltage = @(data, channel, inds) ...
    highpass(double(data.Data.Data(channel,inds))*bitVolts, highPassFreq, fs); % extract voltage from memmapfile, converting to votlage, highpassing, and only return specific channel

% load data
data = memmapfile(file, 'Format', {'int16', [channelNum, smps], 'Data'}, 'Writable', false);


% plot all channels at dft time intervals
figure('color', 'white', 'Units', 'pixels', 'Position', [2000 20 1800 950]);
timeStarts = min(times):timeIntervals:max(times);
traces = nan(channelNum, length(timeStarts), windowSize*fs);

% get traces in parallel
disp('collecting traces...')
parfor i = 1:channelNum
    tracesSub = nan(length(timeStarts), windowSize*fs);
    for j = 1:length(timeStarts)
        startInd = find(times>timeStarts(j));
        inds = startInd:startInd+windowSize*fs-1;
        tracesSub(j,:) = feval(getVoltage, data, i, inds);
    end
    traces(i,:,:) = tracesSub;
end


% plot traces
disp('plotting...')
for i = 1:channelNum
    for j = 1:length(timeStarts)
        subplot(1,length(timeStarts),j); hold on
        plot(timesSub, squeeze(traces(i,j,:)) + find(i==sortInds)*(range(yLims)));
    end
end

% pimp fig
for i = 1:length(timeStarts)
    subplot(1,length(timeStarts),i);
    
    set(gca, 'XLim', [0 timesSub(end)], ...
        'YLim', [0 (channelNum+.5)*range(yLims)], ...
        'Visible', 'off')
    
    text(0, (channelNum+1)*(range(yLims)), ...
        [num2str(round((timeStarts(i)-min(timeStarts))/60)) ' minutes'])
end

% save that ish
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures', 'ephysDrift', [session '.png']));
disp('all done!')


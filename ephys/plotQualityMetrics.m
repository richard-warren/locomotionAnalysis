function plotQualityMetrics(session, cell)

% settings
window = [-.5 1.5]; % ms pre and post spike time to plot
spkNum = 10000;
xSpacing = 4;
ySpacing = 6;
timeBinNum = 4;



% initializations
files = dir(fullfile(getenv('OBSDATADIR'), 'sessions', session));
ephysFolder = files([files.isdir] & contains({files.name}, 'ephys_')).name;
contFiles = dir(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, '*.continuous'));
channelNum = length(contFiles);
fileNameBase = contFiles(1).name(1:3);


% get spike times for good units
addpath(fullfile(getenv('GITDIR'), 'npy-matlab'))
allSpkInds = readNPY(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, 'spike_times.npy'));
clusters = readNPY(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, 'spike_clusters.npy'));
clusterGroups = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, 'cluster_groups.csv'));
unit_ids = clusterGroups.cluster_id(strcmp(clusterGroups.group, 'good'));
spkInds = allSpkInds(clusters==unit_ids(cell));


% get channel mapping
ephysInfo = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', 'ephysInfo.xlsx'), 'Sheet', 'ephysInfo');
mapFile = ephysInfo.map{strcmp(session, ephysInfo.session)};
load(fullfile(getenv('OBSDATADIR'), 'ephys', 'channelMaps', 'kilosort', [mapFile '.mat']), ...
    'xcoords', 'ycoords')


% get fs, microvolts conversion factor, and number of samples
[~, ~, info] = load_open_ephys_data_faster(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, [fileNameBase '_CH1.continuous']));
fs = info.header.sampleRate;
bitVolts = info.header.bitVolts; disp(bitVolts);
file = fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, [fileNameBase '_CHs.dat']);
temp = dir(file);
smps = temp.bytes/2/channelNum; % 2 bytes per sample

% function to extract voltage from binary file
getVoltage = @(data, channel, inds) ...
    double(data.Data.Data(channel,inds))*bitVolts; % extract voltage from memmapfile, converting to votlage, highpassing, and only return specific channel

% load data
data = memmapfile(file, 'Format', {'int16', [channelNum, smps], 'Data'}, 'Writable', false);



% extract waveforms across all channels
spkIndsSub = spkInds(sort(randperm(length(spkInds),min(length(spkInds),spkNum)))); % take random subset of all spikes
spkIndsSub = spkIndsSub(spkIndsSub>range(window)/1000); % ensure spikes aren't too close to the beginning of recording
spkWindow = int64((window(1)/1000*fs) : (window(2)/1000*fs));
spkIndsSubAll = uint64(repmat(spkWindow,1,length(spkIndsSub)) + int64(repelem(spkIndsSub, length(spkWindow)))');

allData = getVoltage(data, 1:channelNum, spkIndsSubAll);
allData = reshape(allData, channelNum, length(spkWindow), []);
allData = permute(allData, [3 1 2]);
allData = double(allData)*bitVolts;
allData = allData - allData(:,:,1); % subtract beginning of trace from rest of trace (a hack of a high pass filter, lol)




% PLOT SPIKE SHAPES ON PROBE
figure('Name', sprintf('%s cell %i', session, cell), 'Color', 'white', 'Position', [2000 50 600 800]); hold on
timeBins = discretize(double(spkIndsSub), timeBinNum);
colors = copper(timeBinNum);

for i = 1:timeBinNum
    for j = 1:channelNum
        trace = squeeze(mean(allData(timeBins==i,j,:),1));
        plot(xcoords(j)*xSpacing + spkWindow, ...
            ycoords(j)*ySpacing + trace, ...
            'Color', colors(i,:), 'LineWidth', 2)
    end
end
set(gca, 'visible', 'off')



% PLOT FIRING RATE + AMPLITUDE OVER TIME
% !!!







% PLOT AMP DISTRIBUTION + AUTOCORRELATION (OVER TIME?)
% !!!





% PLOT FALSE POSITIVES + FALSE NEGATIVES OVER TIME


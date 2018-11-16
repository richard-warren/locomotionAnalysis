% function plotQualityMetrics(session, cell)

% temp
session = '181002_002';
cell = 1;


% settings
window = [-.5 1.5]; % ms pre and post spike time to plot
spkNum = 10000;
xSpacing = 4;
ySpacing = 6;
timeBinNum = 4;



% initializations
meanColor = mean(copper(3),1);
files = dir(fullfile(getenv('OBSDATADIR'), 'sessions', session));
ephysFolder = files([files.isdir] & contains({files.name}, 'ephys_')).name;
contFiles = dir(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, '*.continuous'));
channelNum = length(contFiles);
fileNameBase = contFiles(1).name(1:3);


% get spike times for good units
addpath(fullfile(getenv('GITDIR'), 'npy-matlab'))
addpath(fullfile(getenv('GITDIR'), 'analysis-tools'))
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
[~, timeStamps, info] = load_open_ephys_data_faster(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, [fileNameBase '_CH1.continuous']));
fs = info.header.sampleRate;
bitVolts = info.header.bitVolts;
file = fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, [fileNameBase '_CHs.dat']);
temp = dir(file);
smps = temp.bytes/2/channelNum; % 2 bytes per sample

% function to extract voltage from binary file
getVoltage = @(data, channel, inds) ...
    double(data.Data.Data(channel,inds))*bitVolts; % extract voltage from memmapfile, converting to votlage, highpassing, and only return specific channel

% load data
data = memmapfile(file, 'Format', {'int16', [channelNum, smps], 'Data'}, 'Writable', false);



% extract waveforms across all channels
spkIndsSub = spkInds(round(linspace(2, length(spkInds), spkNum))); % get subpopulation of spikes evenly spaced out
spkWindow = int64((window(1)/1000*fs) : (window(2)/1000*fs));
spkIndsSubAll = uint64(repmat(spkWindow,1,length(spkIndsSub)) + int64(repelem(spkIndsSub, length(spkWindow)))');

allWaveforms = getVoltage(data, 1:channelNum, spkIndsSubAll);
allWaveforms = reshape(allWaveforms, channelNum, length(spkWindow), []);
allWaveforms = permute(allWaveforms, [3 1 2]);
allWaveforms = double(allWaveforms)*bitVolts;
allWaveforms = allWaveforms - allWaveforms(:,:,1); % subtract beginning of trace from rest of trace (a hack of a high pass filter, lol)




%% PLOT SPIKE SHAPES ON PROBE
figure('Name', sprintf('%s cell %i', session, cell), 'Color', 'white', 'Position', [2000 50 600 800]); hold on
timeBins = discretize(double(spkIndsSub), timeBinNum);
colors = copper(timeBinNum);

for j = 1:channelNum
    for i = 1:timeBinNum
        trace = squeeze(mean(allWaveforms(timeBins==i,j,:),1));
        plot(xcoords(j)*xSpacing + spkWindow, ...
            ycoords(j)*ySpacing + trace, ...
            'Color', colors(i,:), 'LineWidth', 2)
    end
    text(double(xcoords(j)*xSpacing+spkWindow(1)), ...
            ycoords(j)*ySpacing, ...
            num2str(j))
end
set(gca, 'visible', 'off')



% PLOT OTHER QUALITY METRICS

% initializations
figure('color', 'white', 'position', [1948 119 1416 823]);
xLims = [0 range(timeStamps)/60];


% PLOT AMP OVER TIME
% settings
circSize = 20;


subplot(3,1,1)

% find best channel
meanWaveform = squeeze(mean(allWaveforms,1));
[~, bestChannel] = max(peak2peak(meanWaveform,2));

amplitudes = peak2peak(squeeze(allWaveforms(:,bestChannel,:)),2);
channelData = getVoltage(data, bestChannel, 1:smps);
chanelData = highpass(channelData, 600, fs);
stdev = std(channelData);

minutes = (timeStamps(spkIndsSub)-timeStamps(1)) / 60;
scatter(minutes, amplitudes/stdev, circSize, copper(length(spkIndsSub)));
set(gca, 'XLim', xLims)
ylabel(sprintf('channel %i SNR', bestChannel))
xlabel('time (min)')


% PLOT FIRING RATE OVER TIME
% settings
timeBins = 1000;
binSmoothing = 10;

subplot(3,1,2)

[counts, edges] = histcounts(timeStamps(spkInds), timeBins);
binWidth = median(diff(edges));
spkRates = counts / binWidth;
spkRates = smooth(spkRates, binSmoothing);
minutes = (edges(1:end-1) + .5*binWidth - timeStamps(1)) / 60; % get bin centers and convert from seconds to minutes
plot(minutes, spkRates, 'LineWidth', 2, 'Color', meanColor)
ylabel('firing rate (Hz)')

set(gca, 'Box', 'off', 'XLim', xLims)



% % AUTOCORRELOGRAM
% 
% % settings
% binWidth = .001;
% width = .04;
% 
% % spksBinary = false(1,length(timeStamps));
% % spksBinary(spkInds) = true;
% spksBinary = histcounts(timeStamps(spkInds), 0:binWidth:timeStamps(end));
% autocorr = xcorr(spksBinary, round(width/2/binWidth), 'coeff');
% autocorr(autocorr==1) = nan;
% 
% % lags = linspace(-width*.5, width*.5, length(autocorr));
% cla; bar(autocorr, 'BarWidth', 1, 'EdgeColor', 'none', 'FaceColor', mean(copper(3),1))


% FALSE POSITIVE RATE
% settings
timeBins = 200;
binSize = 5*60;
refractoryPeriod = .003;
censoredPeriod = .001;

subplot(3,1,3)

fpRates = nan(1,timeBins);
centers = linspace(0, timeStamps(spkInds(end)), timeBins);

for i = 1:timeBins
    
    binMin = max(0, centers(i)-.5*binSize);
    binMax = min(timeStamps(spkInds(end)), centers(i)+.5*binSize);
    
    spkBins = timeStamps(spkInds)>=binMin & timeStamps(spkInds)<binMax;
    violations = sum(diff(timeStamps(spkInds(spkBins)))<refractoryPeriod);
    if violations>0
        rhs = 2 * (refractoryPeriod-censoredPeriod) * length(spkInds)^2 / violations / (binMax-binMin); % what the fuck is this math about??? is this really correct???
        fpRates(i) = .5 - .5*sqrt((rhs-4)/rhs);
    else
        fpRates(i) = 0;
    end
end

plot((centers-timeStamps(1))/60, fpRates, 'LineWidth', 2, 'Color', meanColor)
set(gca, 'box', 'off', 'XLim', xLims, 'YLim', [0 .05])
ylabel('false positive rate')
xlabel('time (min)')









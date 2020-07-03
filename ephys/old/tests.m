
% settings
session = '181019_002';
highPassFreq = 300;
fs = 30000;

% initializations
files = dir(fullfile(getenv('OBSDATADIR'), 'sessions', session));
ephysFolder = files([files.isdir] & contains({files.name}, 'ephys_')).name;
addpath(fullfile(getenv('GITDIR'), 'npy-matlab'))
addpath(fullfile(getenv('GITDIR'), 'analysis-tools'))
npyFiles = dir(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, '*.npy'));
channelNum = length(dir(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, '*.continuous')));

% get fs and microvolts conversion factor
[~, ~, info] = load_open_ephys_data_faster(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, '107_CH1.continuous'));
fileFs = info.header.sampleRate;
bitVolts = info.header.bitVolts;
if fileFs~=fs; disp('sampling freq error! abandon ship!'); keyboard; end

% function to extract voltage from binary file
getVoltage = @(data, channel, inds) ...
    highpass(double(data.Data.Data(channel,inds))*bitVolts, highPassFreq, fs); % extract voltage from memmapfile, converting to votlage, highpassing, and only return specific channel

% get total number of samples
file = fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, '107_CHs.dat');
temp = dir(file);
smps = temp.bytes/2/channelNum; % 2 bytes per sample

% load data
data = memmapfile(file, 'Format', {'int16', [channelNum, smps], 'Data'}, 'Writable', false);



%% plot sample
time = 5; % seconds
close all; figure;
offset = 500;

for i = 1:32
    plot(getVoltage(data,i,1:time*fs) + offset*i); hold on
    text(0, offset*i, num2str(i))
    pause(.01)
end

set(gca, 'Visible', 'off', 'YLim', [0 offset*(i+1)])
pimpFig


%% get binned spike rates over recording duration

% settings
timeToSmp = 60;
channel = 25;
spkThresh = 4;
timeBins = 100;


getSpikeThresh = @(x) -(median(abs(x))/.6745) * spkThresh;
rates = nan(32, timeBins);

for i = 1:32
    disp(i)
    channelTrace = getVoltage(data,channel,1:smps);
    thresh = getSpikeThresh(channelTrace);
    spkBins = logical([0, diff(channelTrace<thresh)==1]);
    times = linspace(0,smps/fs,smps);
    spkTimes = times(spkBins);
    [counts, edges] = histcounts(spkTimes, timeBins);
    rates(i,:) = counts / median(diff(edges));
end


%% sort channels by depth and show rates in heat map

figure;
plot(channelTrace); hold on;
line(get(gca, 'XLim'), [thresh thresh])
scatter(find(spkBins), channelTrace(spkBins))
pimpFig


















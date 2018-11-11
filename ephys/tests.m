
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
getVoltage = @(data, channel) highpass(double(data.Data.Data(channel,:))*info.header.bitVolts, highPassFreq, fs);

% % load all npy files
% npyData = struct();
% for i = 1:length(npyFiles)
%     npyData.(npyFiles(i).name(1:end-4)) = ...
%         readNPY(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, npyFiles(i).name));
% end

% close all; figure
% scatter(npyData.spike_times, npyData.amplitudes)
% figure; histogram(npyData.amplitudes)

%% play with continuous data

% the following puts data into a convenient memory map file 

% create copy of binary file to play with
file = fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, '107_CHs');
copyfile([file '.dat'], [file 'Normalized.dat'])

% get total number of samples
temp = dir([file 'Normalized.dat']);
smps = temp.bytes/2/channelNum; % 2 bytes per sample

% acquire and normalize data
data = memmapfile([file 'Normalized.dat'], ...
    'Format', {'int16', [channelNum, smps], 'Data'}, 'Writable', true);


%% plot all channels for brief window
time = 5; % seconds
close all; figure;
offset = 500;

for i = 1:32
    trace = double(data.Data.Data(i,1:time*fs)) * info.header.bitVolts;
    trace = highpass(double(trace), highPassFreq, fs);
    plot(trace + offset*i); hold on
    text(0, offset*i, num2str(i))
end

set(gca, 'Visible', 'off', 'YLim', [0 offset*(i+1)])
pimpFig


%% plot single trace amplitude metric over time

% settings
% timeToSample = 60; % randomly
channel = 25;
spkThresh = 4;

getSpikeThresh = @(x) -(median(abs(x))/.6745) * spkThresh;
thresh = getSpikeThresh(getVoltage(data,channel));
spkBins = logical([0, diff(getVoltage(data,channel)>thresh)==1]);


figure;
plot(getVoltage(data,channel)); hold on;
line(get(gca, 'XLim'), [thresh thresh])
scatter(find(spkBins), repmat(thresh,1,sum(spkBins)))
pimpFig


















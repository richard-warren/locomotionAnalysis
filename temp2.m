% SHOW TRACES FROM A DATA FILE

% temp
session = '180920_002';
channel = 2;
smpsToShow = [1:(30000*1)] + (10*60*30000);

files = dir(fullfile(getenv('OBSDATADIR'), 'sessions', session));
ephysFolder = files([files.isdir] & contains({files.name}, 'ephys_')).name;
contFiles = dir(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, '*.continuous'));
channelNum = length(contFiles);
fileNameBase = contFiles(1).name(1:3);
file = fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, [fileNameBase '_CHs.dat']);
temp = dir(file);
smps = temp.bytes/2/channelNum; % 2 bytes per sample


% function to extract voltage from binary file
getVoltage = @(data, channel, inds) data.Data.Data(channel,inds); % extract voltage from memmapfile, converting to votlage, highpassing, and only return specific channel

% load data
data = memmapfile(file, 'Format', {'int16', [channelNum, smps], 'Data'}, 'Writable', false);


channelData = getVoltage(data, channel, smpsToShow);

close all; figure('Position', [138 428 1680 411]);
subplot(2,1,1); plot(highpass(double(channelData), 1000, 30000)); ylabel('filtered')
subplot(2,1,2); plot(channelData); ylabel('dat file')

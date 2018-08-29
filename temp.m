fileId = fopen('Z:\RAW\obstacleData\ephys\tests\2018-08-26_15-58-09\100_CHs.dat');
data = fread(fileId, 'int16', 0, 'l');
%%

close all; figure;
% channelsToShow = 1:32;
channelsToShow = [2 14 19 23 28];
fs = 30000;
duration = 1;
colors = copper(length(channelsToShow));

for i = channelsToShow
    subplot(length(channelsToShow),1,find(channelsToShow==i))
    signal = data(i:channelNum:end,:);
    signal = highpass(signal, 300, 30000);
    plot(signal(1:fs*duration), 'color', colors(channelsToShow==i,:));
    set(gca, 'visible', 'off')
end
pimpFig




%%

fid = fopen('Z:\RAW\obstacleData\ephys\tests\2018-08-26_15-58-09\100_CH18.continuous');
hdr = fread(fid, 1024, 'char*1');
timestamp = fread(fid, 1, 'int64',0,'l');
N = fread(fid, 1, 'uint16',0,'l');
recordingNumber = fread(fid, 1, 'uint16', 0, 'l');
samples = fread(fid, N, 'int16',0,'b');
recordmarker = fread(fid, 10, 'char*1');
fclose(fid);
figure; plot(samples);
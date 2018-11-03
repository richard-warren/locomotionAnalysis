

% settings
channels = 1:2:64; % 4
highPass = 300;
fs = 30000;
folder = 'Y:\obstacleData\sessions\181030_000\ephys_2018-10-30_18-27-21';
yLims = [-500 500];
windowSize = .5;
timeIntervals = 5*60;

% mapping
% channelIn = [31 27 22 18 28 23 21 26 29 24 20 25 30 19 32 17 1 16 3 14 9 10 8 2 7 15 11 12 6 13 5 4];
% channelOut = [31 28 25 22 19 16 13 10 7 4 32 26 20 14 8 2 1 5 11 17 23 29 3 6 9 12 15 18 21 24 27 30];
channelIn = 1:64;
channelOut = 1:64;

% initializations
addpath(fullfile(getenv('GITDIR'), 'analysis-tools'))
[~, times] = load_open_ephys_data_faster(fullfile(folder, '107_CH1.continuous'));
timeStarts = min(times):timeIntervals:max(times);

figure('color', 'white', 'Units', 'normalized', 'Position', [0 0 1 1]);

for i = 1:length(channels)
    disp(i)
    chan = channelOut(find(channelIn==i,1,'first'));
    
    [data, times, info] = load_open_ephys_data_faster(fullfile(folder, ['107_CH' num2str(channels(i)) '.continuous']));
    dataHpf = highpass(data, highPass, fs);
    
    for j = 1:length(timeStarts)
        bins = times>timeStarts(j) & times<timeStarts(j)+windowSize;
        subplot(1,length(timeStarts),j);
        plot(times(bins), dataHpf(bins) + (32-chan-1)*(range(yLims))); hold on;
    end
    pause(.5)
end


% pimp fig
for i = 1:length(timeStarts)
    subplot(1,length(timeStarts),i);
    bins = times>timeStarts(i) & times<timeStarts(i)+windowSize;
    set(gca, 'XLim', [min(times(bins)) max(times(bins))], ...
        'YLim', [yLims(1) (yLims(1)+range(yLims)*length(channels))], ...
        'Visible', 'off')
    text(min(times(bins)), 32*(range(yLims)), [num2str(round((timeStarts(i)-min(timeStarts))/60)) ' minutes'])
end

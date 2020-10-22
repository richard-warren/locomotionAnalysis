function makeTouchSensorVid2(session, trialNum)

% this version is modified for a single touch sensor, not the two access cantelever force sensor

% settings
editedDir = [getenv('OBSDATADIR') 'editedVid\'];
fps = 25;
sensorXWidth = .5; % seconds
yLims = [-.5 1];
vidScaling = 1.5;

% load session data
vidTop = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);
load([getenv('OBSDATADIR') 'sessions\' session '\run.mat'], 'touch');
load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], 'obsOnTimes', 'obsOffTimes', 'frameTimeStamps', 'obsPixPositions');

% initializations
trials = randperm(length(obsOnTimes), min(trialNum, length(obsOnTimes)));
trials = sort(trials);
% trials = 30:40; % temp
vidWriter = VideoWriter(sprintf('%ssensorVid%s', editedDir, session), 'MPEG-4');
set(vidWriter, 'FrameRate', fps)
open(vidWriter)


% prepare figure
fig = figure('position', [2000 0 vidTop.Width*vidScaling vidTop.Height*2*vidScaling], 'menubar', 'none');

subplot(2,1,1);
framePreview = imshow(read(vidTop,1));
set(gca, 'position', [0 .5 1 .5])

subplot(2,1,2);
trace = plot(1, 1, 'color', 'white', 'linewidth', 1); hold on;
set(gca, 'ylim', yLims, 'position', [0 0 1 .5], 'box', 'off', 'xtick', [], 'ytick', [], 'color', 'black')

w = waitbar(0, 'editing video...');

for i = 1:length(obsOnTimes)%trials
    
    % get trial inds
    trialInds = find(frameTimeStamps>=obsOnTimes(i) & frameTimeStamps<=obsOffTimes(i) & ...
        obsPixPositions'>0 & obsPixPositions'<=vidTop.Width);
    % damn son, the following line is a major hack:
    trialInds = cat(1, [trialInds(1)-50:trialInds(1)-1]', trialInds);
    
    
    for j = trialInds'
        
        % get frame
        frame = read(vidTop, j);
        frame = insertText(frame, [size(frame,2) size(frame,1)], num2str(i),...
                           'BoxColor', 'black', 'AnchorPoint', 'RightBottom', 'TextColor', 'white');
        
        % get sensor traces
        sensorBins = touch.times>frameTimeStamps(j)-sensorXWidth & touch.times<=frameTimeStamps(j);
        top = touch.values(sensorBins);
        times = linspace(-sensorXWidth,0,length(top));
        
        % update frame and plot
        set(framePreview, 'CData', frame);
        set(trace, 'XData', times, 'YData', top);
        set(gca, 'xlim', [times(1), times(end)])
        
        % write to video
        writeVideo(vidWriter, getframe(gcf));
        
    end
    
    % update waitbar
    waitbar(i/length(obsOnTimes))
    
end

close(vidWriter)
close(fig)
close(w)

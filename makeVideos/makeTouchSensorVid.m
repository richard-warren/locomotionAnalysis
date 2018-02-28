function makeTouchSensorVid(session, trials)

% settings
editedDir = [getenv('OBSDATADIR') 'editedVid\'];
fps = 15;
sensorXWidth = .5; % seconds
yLims = [-.3 .3];
vidScaling = 1.5;
smoothSmps = 1;

% load session data
vidTop = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);
load([getenv('OBSDATADIR') 'sessions\' session '\run.mat'], 'backFron', 'topDown');
load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], 'obsOnTimes', 'obsOffTimes', 'frameTimeStamps', 'obsPixPositions');

% initializations
trials = randperm(length(obsOnTimes), min(trials, length(obsOnTimes)));
trials = sort(trials);
trials = 30:40; % temp
vidWriter = VideoWriter(sprintf('%ssensorVid%s', editedDir, session), 'MPEG-4');
set(vidWriter, 'FrameRate', fps)
open(vidWriter)

% normalize and smooth sensor values
backFron.values = smooth(backFron.values - mean(backFron.values), smoothSmps);
topDown.values = smooth(topDown.values - mean(topDown.values), smoothSmps);

% prepare figure
fig = figure('position', [2000 0 vidTop.Width*vidScaling vidTop.Height*2*vidScaling], 'menubar', 'none');

subplot(2,1,1);
framePreview = imshow(read(vidTop,1));
set(gca, 'position', [0 .5 1 .5])

subplot(2,1,2);
topTrace = plot(1, 1, 'color', [51, 204, 255] / 255, 'linewidth', 3); hold on;
frontTrace = plot(1, 1, 'color', [102, 255, 102] / 255, 'linewidth', 3);
set(gca, 'ylim', yLims, 'position', [0 0 1 .5], 'box', 'off', 'xtick', [], 'ytick', [], 'color', 'black')
legend({'up/down', 'front/back'}, 'location', 'northwest', 'textcolor', 'white')
legend('boxoff')

w = waitbar(0, 'editing video...');

for i = trials
    
    % get trial inds
    trialInds = find(frameTimeStamps>=obsOnTimes(i) & frameTimeStamps<=obsOffTimes(i) & ...
        obsPixPositions'>0 & obsPixPositions'<=vidTop.Width);
    
    for j = trialInds'
        
        % get frame
        frame = read(vidTop, j);
        frame = insertText(frame, [size(frame,2) size(frame,1)], num2str(i),...
                           'BoxColor', 'black', 'AnchorPoint', 'RightBottom', 'TextColor', 'white');
        
        % get sensor traces
        sensorBins = topDown.times>frameTimeStamps(j)-sensorXWidth & topDown.times<=frameTimeStamps(j);
        top = topDown.values(sensorBins);
        front = backFron.values(sensorBins);
        times = linspace(-sensorXWidth,0,length(top));
        
        % update frame and plot
        set(framePreview, 'CData', frame);
        set(topTrace, 'XData', times, 'YData', top);
        set(frontTrace, 'XData', times, 'YData', front);
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

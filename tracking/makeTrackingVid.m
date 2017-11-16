function makeTrackingVid

% makes video of top and bot views with tracked positions overlayed
% if duration==0 the video is made for the whole file

% user settings
fps = 25; % playback fps
circSize = 75;

% load data
load('C:\Users\rick\Google Drive\columbia\obstacleData\svm\trackedData\locationsBot.mat', 'locationsBot')
load('C:\Users\rick\Google Drive\columbia\obstacleData\svm\trackedData\locationsTop.mat', 'locationsTop')
% locationsBot = fixTracking(locationsBot);
% locationsTop = fixTracking(locationsTop);
% locationsTop.x = locationsBot.x; % replace top x values with those from bottom, which are more reliable

vidTop = VideoReader('C:\Users\rick\Google Drive\columbia\obstacleData\svm\testVideo\runTop.mp4');
vidBot = VideoReader('C:\Users\rick\Google Drive\columbia\obstacleData\svm\testVideo\runBot.mp4');
frameTop = read(vidTop,1);
frameBot = read(vidBot,1);


% initializations
vidWrite = VideoWriter('C:\Users\rick\Google Drive\columbia\obstacleData\svm\testVideo\trackingSample', 'MPEG-4');
set(vidWrite, 'FrameRate', fps);
open(vidWrite);
colors = winter(4);


% set up figure
close all
mainFig = figure('color', [0 0 0], 'position', [1925, 50, vidTop.width, (vidBot.Height + vidTop.Height)], 'menubar', 'none');

topAxis = subplot(2,1,1, 'units', 'pixels');
colormap gray
topImshow = image(frameTop, 'parent', topAxis, 'CDataMapping', 'scaled'); hold on;
set(topAxis, 'visible', 'off', 'CLim', [0 255])
scatterTop = scatter(topAxis, locationsTop.x(1,:), locationsTop.z(1,:), circSize, colors, 'filled');
txt = text(topAxis, 5 , 5, '', 'color', [1 1 1]);

botAxis = subplot(2,1,2, 'units', 'pixels');
botImshow = image(frameBot, 'parent', botAxis, 'CDataMapping', 'scaled'); hold on;
scatterBot = scatter(botAxis, locationsBot.x(1,:), locationsBot.y(1,:), circSize, colors, 'filled');
set(botAxis, 'visible', 'off', 'CLim', [0 255])


set(topAxis, 'position', [0 vidBot.Height vidTop.Width vidTop.Height]);
set(botAxis, 'position', [0 0 vidBot.Width vidBot.Height]);


% write video
for i = 1:vidTop.NumberOfFrames
    
    % plot top frame
    frameTop = read(vidTop,i);
    frameTop = imadjust(squeeze(frameTop(:,:,1)), [.1 1], [0 1]);
    set(topImshow, 'CData', frameTop);
    set(txt, 'String', ['f = ' num2str(i)])
%     keyboard

    % plot bot frame
    frameBot = read(vidBot,i);
    frameBot = imadjust(squeeze(frameBot(:,:,1)), [.1 1], [0 1]);
    frameBot = medfilt2(frameBot, [4 8], 'symmetric');
    set(botImshow, 'CData', frameBot);

    % update circles
    set(scatterTop, 'XData', locationsTop.x(i,:), 'YData', locationsTop.z(i,:))
    set(scatterBot, 'XData', locationsBot.x(i,:), 'YData', locationsBot.y(i,:))

    writeVideo(vidWrite, getframe(gcf));
end

close(vidWrite);
close(mainFig);



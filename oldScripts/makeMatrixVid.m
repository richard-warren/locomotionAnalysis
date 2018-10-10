% function makeMatrixVid
% make sweet vid with mouse running at real splayback speed then going super slow mo at moment of wisk contact, and drawing kinematic traces on top of image // dope

% settings
session = '180124_001';
trials = [4 5 10 11 14 23 41 43 63 71]; % 5 14 43
prePostTime = [-.2 .2]; % (s) time to add to beginning and end of a trial (before and after obs is engaged) // seconds are in real time
slowSpeed = .05; % how much to slow down vid after wisk contacts obstacle
contrastLims = [.2 1]; % pixels at these proportional values are mapped to 0 and 255
cmap = 'hsv';
circSize = 150;
circs = 4;
circSeparation = 2; % how many circs separate each frame
trailDarkening = .5;
iti = 0; % time to put black frames in between trials

% initializations
circIndOffsets = -(circs-1)*circSeparation : circSeparation : 0;
colors = eval([cmap '(4)']);
colors = repelem(colors,circs,1);
trailBins = logical(mod((1:4*circs),circs)~=0);
colors(trailBins,:) = colors(trailBins,:) * trailDarkening;
circSizes = [linspace(circSize*.1, circSize*.5, circs-1) circSize];
circSizes = repmat(circSizes,1,4);

load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], ...
    'obsOnTimes', 'obsOffTimes', 'frameTimeStamps');
load([getenv('OBSDATADIR') 'sessions\' session '\tracking\locationsBotCorrected.mat'], 'locations');
load([getenv('OBSDATADIR') 'sessions\' session '\wiskContactTimes.mat'], 'contactTimes');
vidTop = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);
vidBot = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4']);
vidWriter = VideoWriter([getenv('OBSDATADIR') 'editedVid\' sprintf('matrixVid%sspeed&%.2ftrials%s', session, slowSpeed, num2str(trials))], 'MPEG-4');
frameRate = 25; % frameRate of video file to be written
set(vidWriter, 'FrameRate', frameRate);
open(vidWriter);

% set up figure
mainFig = figure('name', session, 'color', [0 0 0], 'position', [1925, 50, vidTop.width, (vidBot.Height + vidTop.Height)], 'menubar', 'none');
colormap gray

topAxis = subplot(2,1,1, 'units', 'pixels');
topIm = image(rgb2gray(read(vidTop,1)), 'parent', topAxis, 'CDataMapping', 'scaled'); hold on;
set(topAxis, 'visible', 'off', 'CLim', [0 255])
scatterTop = scatter(topAxis, zeros(1,4*circs), zeros(1,4*circs), circSizes, colors, 'filled');

botAxis = subplot(2,1,2, 'units', 'pixels');
botIm = image(rgb2gray(read(vidBot,1)), 'parent', botAxis, 'CDataMapping', 'scaled'); hold on;
scatterBot = scatter(botAxis, zeros(1,4*circs), zeros(1,4*circs), circSizes, colors, 'filled');
set(botAxis, 'visible', 'off', 'CLim', [0 255])

set(topAxis, 'position', [0 vidBot.Height vidTop.Width vidTop.Height]);
set(botAxis, 'position', [0 0 vidBot.Width vidBot.Height]);



% iterate through trials
for i = 1:length(trials)
    
    % iterate through fast frames (prior to wisk contact)
    fastBins = frameTimeStamps>=(obsOnTimes(trials(i))+prePostTime(1)) & ...
                    frameTimeStamps<=contactTimes(trials(i));
    fastInds = find(fastBins);
    trialFrameTimes = frameTimeStamps(fastInds(1)) : (1/frameRate) : frameTimeStamps(fastInds(end)); % these are the 'desired' frame times given the specified frameRate -- will find the frames closest to the desired frames
    trialFastInds = fastInds(knnsearch(frameTimeStamps(fastInds), trialFrameTimes'));
    
    slowBins = frameTimeStamps>contactTimes(trials(i)) & ...
                frameTimeStamps<=(obsOffTimes(trials(i))+prePostTime(2));
    slowInds = find(slowBins);
    trialFrameTimes = frameTimeStamps(slowInds(1)) : (slowSpeed/frameRate) : frameTimeStamps(slowInds(end)); % these are the 'desired' frame times given the specified frameRate -- will find the frames closest to the desired frames
    trialSlowInds = slowInds(knnsearch(frameTimeStamps(slowInds), trialFrameTimes'));
    
    trialInds = [trialFastInds; trialSlowInds];
    
    
    
    for j = 1:length(trialInds)
        
        % get top frame
        frameTop = rgb2gray(read(vidTop, trialInds(j)));
        frameTop = imadjust(frameTop, contrastLims, [0 1]);
        set(topIm, 'CData', frameTop);
        
        % get bot frame
        frameBot = rgb2gray(read(vidBot, trialInds(j)));
        frameBot = imadjust(frameBot, contrastLims, [0 1]);
        set(botIm, 'CData', frameBot);
        
        % add scatter points
        x = squeeze(locations.locationsCorrected(trialInds(j)+circIndOffsets,1,:));
        y = squeeze(locations.locationsCorrected(trialInds(j)+circIndOffsets,2,:));
        set(scatterBot, 'XData', x(:), 'YData', y(:));
        
        % write to video
        frame = getframe(gcf);
        writeVideo(vidWriter, frame);
    end
    
    % add black frames between trials
    for j = 1:(iti*frameRate)
        writeVideo(vidWriter, zeros(size(frame.cdata)));
    end
end


close(vidWriter)
close(mainFig)
disp('all done')
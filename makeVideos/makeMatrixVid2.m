function makeMatrixVid2(filename, session, trials, allSlowMotion, slowSpeed, dropEveryOtherFrame)
% make sweet vid with mouse running at real splayback speed then going super slow mo at moment of wisk contact, and drawing kinematic traces on top of image // dope
% if allSlowMotion true, then entire video is played back at slow speed

% settings
% prePostTime = [-.2 .2]; % (s) time to add to beginning and end of a trial (before and after obs is engaged) // seconds are in real time
prePostTime = [-.1 .1]; % (s) time to add to beginning and end of a trial (before and after obs is engaged) // seconds are in real time
contrastLims = [.2 1]; % pixels at these proportional values are mapped to 0 and 255
cmap = 'hsv';
circSize = 60;
circs = 4;
circSeparation = 2; % how many circs separate each frame
trailDarkening = .5;
iti = 0; % time to put black frames in between trials
featuresToOmmit = {'obsHigh_bot', 'obsLow_bot', 'obs_top'};
if ~exist('dropEveryOtherFrame', 'var'); dropEveryOtherFrame=false; end



% initializations


% load data
load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], ...
    'obsOnTimes', 'obsOffTimes', 'frameTimeStamps', 'wiskContactTimes');
locationsTable = readtable([getenv('OBSDATADIR') 'sessions\' session '\trackedFeaturesRaw.csv']); % get raw tracking data
[locations, features] = fixTrackingDLC(locationsTable, frameTimeStamps);
featuresToShowBins = ~ismember(features, featuresToOmmit);
features = features(featuresToShowBins);

% get circle size vector
circIndOffsets = -(circs-1)*circSeparation : circSeparation : 0;
circSizes = [linspace(circSize*.1, circSize*.5, circs-1) circSize];
circSizes = repmat(circSizes,1,length(features));

% get colors s.t. same feature in top and bot is same color
[uniqueFeaturs, ~, inds] = unique(cellfun(@(x) x(1:end-4), features, 'UniformOutput', false)); % list of features that appear in top and bot views
colors = eval([cmap '(' num2str(length(uniqueFeaturs)) ')']);
colors = colors(randperm(size(colors,1)),:); % this is a hack - randomizes color order so that colors of paws are not all similar to one another
colors = colors(inds, :);
colors = repelem(colors,circs,1);
trailBins = logical(mod((1:length(features)*circs),circs)~=0);
colors(trailBins,:) = colors(trailBins,:) * trailDarkening;


vidTop = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);
vidBot = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4']);
vidWriter = VideoWriter(filename, 'MPEG-4');
frameRate = 50; % frameRate of video file to be written
if dropEveryOtherFrame; frameRate = frameRate/2; end
set(vidWriter, 'FrameRate', frameRate);
open(vidWriter);

% set up figure
fig = figure('name', session, 'color', [0 0 0], 'position', [1925, 50, vidTop.width, (vidBot.Height + vidTop.Height)], 'menubar', 'none');
ax = axes('position', [0 0 1 1], 'CLim', [0 255]);
colormap gray

frame = cat(1, rgb2gray(read(vidTop,1)), rgb2gray(read(vidBot,1)));
im = image(frame, 'CDataMapping', 'scaled'); hold on;
scat = scatter(zeros(1,length(features)*circs), zeros(1,length(features)*circs), circSizes, colors, 'filled');
set(ax, 'visible', 'off')



% iterate through trials
for i = 1:length(trials)
    
    % iterate through fast frames (prior to wisk contact)
    fastInds = find(frameTimeStamps>=(obsOnTimes(trials(i))+prePostTime(1)) & ...
                    frameTimeStamps<=wiskContactTimes(trials(i)));
    if dropEveryOtherFrame; fastInds = fastInds(1:2:length(fastInds)); end
    
    trialFrameTimes = frameTimeStamps(fastInds(1)) : (1/frameRate) : frameTimeStamps(fastInds(end)); % these are the 'desired' frame times given the specified frameRate -- will find the frames closest to the desired frames
    trialFastInds = fastInds(knnsearch(frameTimeStamps(fastInds), trialFrameTimes'));
    
    if allSlowMotion
        slowInds = find(frameTimeStamps>=(obsOnTimes(trials(i))+prePostTime(1)) & ...
                   frameTimeStamps<=(obsOffTimes(trials(i))+prePostTime(2)));
        trialFastInds = [];
    else
        slowInds = find(frameTimeStamps>wiskContactTimes(trials(i)) & ...
                frameTimeStamps<=(obsOffTimes(trials(i))+prePostTime(2)));
        if dropEveryOtherFrame
            slowInds = slowInds(1:2:length(slowInds));
        end
    end
    trialFrameTimes = frameTimeStamps(slowInds(1)) : (slowSpeed/frameRate) : frameTimeStamps(slowInds(end)); % these are the 'desired' frame times given the specified frameRate -- will find the frames closest to the desired frames
    trialSlowInds = slowInds(knnsearch(frameTimeStamps(slowInds), trialFrameTimes'));
    
    trialInds = [trialFastInds; trialSlowInds];
    
    
    
    for j = 1:length(trialInds)
        
        % get frame
        frame = cat(1, rgb2gray(read(vidTop,trialInds(j))), rgb2gray(read(vidBot,trialInds(j))));
        frame = imadjust(frame, contrastLims, [0 1]);
        set(im, 'CData', frame);
        
        % add scatter points
        x = squeeze(locations(trialInds(j)+circIndOffsets,1,featuresToShowBins));
        y = squeeze(locations(trialInds(j)+circIndOffsets,2,featuresToShowBins));
        set(scat, 'XData', x(:), 'YData', y(:));
        
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
close(fig)
disp('all done!')
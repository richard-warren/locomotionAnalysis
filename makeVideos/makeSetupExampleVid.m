function makeSetupExampleVid(filename, session, trials)
% this make a wide vid with top view plus drawing of obs off screen, intended to explain how obs appears and is controlled by the speed of the wheel
% should use trials in which there is variability in wheel speed, and ideally some in which mouse moves both forwards and backwards

% settings
prePostTime = [-.1 0]; % (s) time to add to beginning and end of a trial (before and after obs is engaged)
playbackSpeed = 0.15;
fps = 250;
obsRadius = 8;
constrastLims = [.2 1]; % pixels at these proportional values are mapped to 0 and 255
obsFadePixels = 100; % when obs enter right side of vid frame, fade out drawing of obs over the course of this many pixels
obsXOffset = 12; % draw obs this many pixels to the right to account forthe pix positions being based on the bottom rather than the top view

% initializations
load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], ...
    'obsPixPositions', 'obsOnTimes', 'obsOffTimes', 'frameTimeStamps');
locationsTable = readtable([getenv('OBSDATADIR') 'sessions\' session '\trackedFeaturesRaw.csv']);
obsTopY = locationsTable.obs_top_1;
clear locationsTable
vid = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);
frameDims = [vid.Height round(max(obsPixPositions)*.9)];


vidWriter = VideoWriter([getenv('OBSDATADIR') 'editedVid\' filename '.mp4'], 'MPEG-4');
set(vidWriter, 'FrameRate', round(fps*playbackSpeed));
open(vidWriter);

% define lookup vector for obstacle opacity (it will fade out as it enters right side of screen)
transparencyGradient = ones(1, frameDims(2)*2); % make this unnecessarily long, trailing ones at the end
gradInds = round(vid.Width-.5*obsFadePixels) : round(vid.Width+.5*obsFadePixels);
transparencyGradient(1:gradInds(1)) = 0;
transparencyGradient(gradInds) = linspace(0,1,length(gradInds));


% iterate through trials
for i = 1:length(trials)
    disp(i/length(trials))
    
    trialBins = frameTimeStamps>=(obsOnTimes(trials(i))+prePostTime(1)) & ...
                frameTimeStamps<=(obsOffTimes(trials(i))+prePostTime(2));
    obsHeight = median(obsTopY(trialBins' & obsPixPositions>0 & obsPixPositions<vid.Width));
    trialInds = find(trialBins);
    
    % iterate through frames within trial
    for j = 1:length(trialInds)
        
        % get vid frame and incorporate into wide frame
        vidFrame = rgb2gray(read(vid, trialInds(j)));
        vidFrame = imadjust(vidFrame, constrastLims, [0 1]);
        frame = uint8(zeros(frameDims));
        frame(:,1:vid.Width) = vidFrame;
        
        % add obs to frame
        pixPos = round(obsPixPositions(trialInds(j)))+obsXOffset;
        if pixPos>0 && pixPos<(frameDims(2)+2*obsRadius)
            opacity = transparencyGradient(pixPos);
            frame = insertShape(frame, 'FilledCircle', ...
                [pixPos obsHeight+obsRadius obsRadius], ...
                'opacity', opacity);
        end
        
        
        % write to video
        writeVideo(vidWriter, frame);
    end
end


close(vidWriter)
disp('all done')
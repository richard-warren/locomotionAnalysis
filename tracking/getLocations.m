function getLocations(session, makeVid)

% performs paw tracking for a single session!
% if makeVid is true, it also automatically makes a video of the tracking

% settings
% session = '180122_000';
minVel = .6; % use .4
obsPrePost = [.2 .2]; % include this many meters before and after obs turns on
velPositions = [-.08 .08] + 0.3820; % compute trials velocity between these obstacle positions
anchorPtsBot = {[0 0], [1 0], [1 1], [0 1]}; % LH, LF, RF, RH // each entry is x,y pair measured from top left corner
colors = hsv(4); % red green blue purple




% initializations
trackingDir = [getenv('OBSDATADIR') 'sessions\' session '\tracking'];
if ~exist(trackingDir, 'dir'); mkdir(trackingDir); end % make tracking directory if it doesn't already exist

xMapping = [getenv('GITDIR') 'locomotionAnalysis\xAlignment\xLinearMapping.mat'];
load(xMapping, 'xLinearMapping');

load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], 'obsPixPositions', 'frameTimeStamps',...
    'wheelPositions', 'wheelTimes', 'obsPositions', 'obsTimes', 'obsOnTimes', 'obsOffTimes', 'mToPixMapping', 'targetFs')
obsPositions = fixObsPositions(obsPositions, obsTimes, obsOnTimes);

[frameInds, trialVels] = getTrialFrameInds(minVel, obsPrePost, velPositions, frameTimeStamps,...
    wheelPositions, wheelTimes, obsPositions, obsTimes, obsOnTimes, obsOffTimes);

mToPixFactor = median(mToPixMapping(:,1)); % get mapping from meters to pixels

vidBot = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4']);
vidTop = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);



%% get bot potential locations


% settings
scoreThresh = -1;
showTracking = 0;
model1 = [getenv('OBSDATADIR') 'svm\classifiers\pawBot1'];

% svm1
load(model1, 'model', 'subFrameSize');
model1 = model; clear model; subFrameSize1 = subFrameSize; clear subFrameSize;
load([getenv('OBSDATADIR') 'svm\classifiers\pawBot2AlexNet.mat'], 'convNetwork', 'subFrameSize');
model2 = convNetwork; clear convNetwork; subFrameSize2 = subFrameSize; clear subFrameSize;
classNum = model2.Layers(end).OutputSize - 1; % not including NOT PAW class

tic
potentialLocationsBot = getPotentialLocationsBot(vidBot, model1, model2, classNum, ...
    subFrameSize1, subFrameSize2, scoreThresh, obsPixPositions, frameInds, showTracking);
                                                  
save([getenv('OBSDATADIR') 'sessions\' session '\tracking\potentialLocationsBot.mat'], 'potentialLocationsBot');
fprintf('%s: potentialLocationsBot analyzed in %.1f minutes\n', session, (toc/60))


%% get bot locations

% settings
showPotentialLocations = 0;

% initializations
% vidBot = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4']);

locationsBot = getLocationsBot(potentialLocationsBot, anchorPtsBot, frameTimeStamps, vidBot.Width, vidBot.Height, frameInds);
% showLocations(vidBot, frameInds, potentialLocationsBot, fixTracking(locationsBot), showPotentialLocations, .02, anchorPtsBot, colors);
save([getenv('OBSDATADIR') 'sessions\' session '\tracking\locationsBot.mat'], 'locationsBot');
fprintf('%s: locationsBot analyzed\n', session)


%% get potential locations for top (svm)

% settings
scoreThresh = 0;
showTracking = 0;

% initializations
load([getenv('OBSDATADIR') 'svm\classifiers\pawTop1'], 'model', 'subFrameSize');
model1 = model; clear model; subFrameSize1 = subFrameSize; clear subFrameSize;
load([getenv('OBSDATADIR') 'svm\classifiers\pawTop2AlexNet'], 'convNetwork', 'subFrameSize')
model2 = convNetwork; clear convNetwork; subFrameSize2 = subFrameSize;
classNum = model2.Layers(end).OutputSize - 1; % not including NOT PAW class


tic; potentialLocationsTop = getPotentialLocationsTop(vidTop, locationsBot, model1, model2, ...
    classNum, subFrameSize1, subFrameSize2, scoreThresh, frameInds, 1:4, showTracking);
save([getenv('OBSDATADIR') 'sessions\' session '\tracking\potentialLocationsTop.mat'], 'potentialLocationsTop');
fprintf('%s: potentialLocationsTop analyzed in %.2f minutes\n', session, toc/60)

%% get locations for top

% settings
showPotentialLocations = 0;
fs = 250;

% fix x alignment for bottom view
locationsBotFixed = fixTracking(locationsBot);
locationsBotFixed.x = locationsBotFixed.x*xLinearMapping(1) + xLinearMapping(2);


locationsTop = getLocationsTop(potentialLocationsTop, locationsBotFixed,...
    frameInds, wheelPositions, wheelTimes, targetFs, mToPixFactor, obsPixPositions, frameTimeStamps, 1:4, fs);
% showLocations(vidTop, frameInds, potentialLocationsTop, fixTracking(locationsTop),...
%     showPotentialLocations, .02, anchorPtsBot, colors, locationsBotFixed);
save([getenv('OBSDATADIR') 'sessions\' session '\tracking\locationsTop.mat'], 'locationsTop');
fprintf('%s: locationsTop analyzed\n', session)


%% make tracking vid

if makeVid; makeTrackingVid(session, frameInds); end





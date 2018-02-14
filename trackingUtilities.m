
% perform paw tracking for multiple sessions

sessionDirs = uigetdir2([getenv('OBSDATADIR') 'sessions\'], 'select folders to analyze');
%%
steps = {'stance'};
minVel = .4;

for i = 1:length(sessionDirs)
    
    % get locations
    nameInd = find(sessionDirs{i}=='\',1,'last');
%     try
        getLocations(sessionDirs{i}(nameInd+1:end), steps, minVel, false);
%     catch
%         fprintf('FAILED TO ANALYZE %s\n', sessionDirs{i}(nameInd+1:end));
%     end
    
end

%% getLocations on single session

session = '180122_000';
steps = {'top'};

showTracking = true;
minVel = .4; % minVel only applies to potentialLocationsBot analysis // subsequent stages analyze whatever has already been analyzed in potentialLocationsBot
getLocations(session, steps, minVel, showTracking);

%% show tracking for session

% settings
session = '180122_000';
view = 'Top';
showCorrected = 0;
frameDelay = .04;

load([getenv('OBSDATADIR') 'sessions\' session '\tracking\potentialLocations' view '.mat'])
if showCorrected
    load([getenv('OBSDATADIR') 'sessions\' session '\tracking\locations' view 'Corrected.mat'])
    locationsTemp = locations.locationsCorrected;
else
    load([getenv('OBSDATADIR') 'sessions\' session '\tracking\locations' view '.mat'])
%     locationsTemp = fixTracking(locations.locationsRaw, locations.trialIdentities);
    locationsTemp = locations.locationsRaw;
end

vid = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\run' view '.mp4']);
anchorPtsBot = {[0 0], [1 0], [1 1], [0 1]}; % LH, LF, RF, RH // each entry is x,y pair measured from top left corner


showLocations(vid, find(locations.isAnalyzed), eval(['potentialLocations' view]), ...
    locationsTemp, locations.trialIdentities, 0, frameDelay, anchorPtsBot, hsv(4));

%% exclude trials

session = '180125_001';
trials = [54 63 65 73 74 88 91 95 97 98 107];

setTrialExclusion(session, trials);


%% correct tracking

% settings
session = '180123_003';
view = 'Bot';
obsCenter = .382;
obsPrePostTop = [.2 .1];
minVel = .4; velPrePost = [.08 .08];

outputFile = [getenv('OBSDATADIR') 'sessions\' session '\tracking\locations' view 'Corrected.mat'];
if exist(outputFile, 'file')
    load(outputFile, 'locations')
else
    load([getenv('OBSDATADIR') 'sessions\' session '\tracking\locations' view '.mat'], 'locations')
end

vid = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\run' view '.mp4']);
frameInds = find(locations.isAnalyzed);
anchorPts = {[0 0], [1 0], [1 1], [0 1]}; % LH, LF, RF, RH

switch view
    case 'Bot'
        frameDelay = .025;
        correctTracking(outputFile, vid, locations, frameInds, frameDelay, anchorPts);
    case 'Top'
        frameDelay = .1;
        paws = [2 3];
        load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], 'obsPixPositions', 'frameTimeStamps',...
            'wheelPositions', 'wheelTimes', 'obsPositions', 'obsTimes', 'obsOnTimes', 'obsOffTimes')
        obsPositions = fixObsPositions(obsPositions, obsTimes, obsOnTimes);
        frameIndsTop = getTrialFrameInds(minVel, obsCenter, obsPrePostTop, velPrePost, frameTimeStamps,...
            wheelPositions, wheelTimes, obsPositions, obsTimes, obsOnTimes, obsOffTimes);

        correctTracking(outputFile, vid, locations, frameIndsTop, frameDelay, anchorPts, paws);
end

%% get avg trial stats

sessionDirs = uigetdir2([getenv('OBSDATADIR') 'sessions\'], 'select folders to analyze');
avgFramesPerTrial = nan(1,length(sessionDirs));

for i = 1:length(sessionDirs)
    
    % get locations
    nameInd = find(sessionDirs{i}=='\',1,'last');
    session = sessionDirs{i}(nameInd+1:end);
    load([getenv('OBSDATADIR') 'sessions\' session '\tracking\velocityInfo.mat'], 'trialVels')
    load([getenv('OBSDATADIR') 'sessions\' session '\tracking\locationsBot.mat'], 'locations')
    avgFramesPerTrial(i) = sum(locations.isAnalyzed) / length(trialVels);
    
end


%% view pose regression training examples



% initializations
imgDir = [getenv('TRAININGEXAMPLESDIR') 'poseRegression\fullRes\'];
originalImSize = [230 396];
load([imgDir 'pawLocations.mat'], 'features')

close all; figure();
img = imread([getenv('TRAININGEXAMPLESDIR') 'poseRegression\fullRes\imgs\img1.tif']);
preview = imshow(img);
hold on; scatterTruth = scatter(gca, [0 0 0 0], [0 0 0 0], 100, hsv(4), 'filled');
% hold on; scatterPredict = scatter(gca, [0 0 0 0], [0 0 0 0], 200, hsv(4));
set(gcf, 'position', [646   173   984   750]);

% pimpFig;


for imNum = randperm(height(features), 20)
    disp(imNum)

    % get image and make predictions!
    img = imread([getenv('TRAININGEXAMPLESDIR') 'poseRegression\fullRes\imgs\img' num2str(imNum) '.tif']);

    % show results
    set(preview, 'CData', img)
    set(scatterTruth, 'XData', table2array(features(imNum, [1 3 5 7]+1))*size(img,2), ...
        'YData', table2array(features(imNum, [2 4 6 8]+1))*size(img,1));
    
    pause(2)
end

%% get getFrameInds

session = '180122_000';

obsCenter = .382; % where obs is at center of wheel
obsPrePostBot = [.55 .25];
obsPrePostTop = [.2 .12];
velPrePost = [.08 .08]; % compute trials velocity between these obstacle positions
anchorPtsBot = {[0 0], [1 0], [1 1], [0 1]}; % LH, LF, RF, RH // each entry is x,y pair measured from top left corner

load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], 'obsPixPositions', 'frameTimeStamps',...
    'wheelPositions', 'wheelTimes', 'obsPositions', 'obsTimes', 'obsOnTimes', 'obsOffTimes', 'mToPixMapping', 'targetFs')
load([getenv('OBSDATADIR') 'sessions\' session '\tracking\potentialLocationsBot.mat'], 'potentialLocationsBot')
load([getenv('OBSDATADIR') 'sessions\' session '\tracking\locationsBotCorrected.mat'], 'locations')
obsPositions = fixObsPositions(obsPositions, obsTimes, obsOnTimes);

[frameIndsBot, trialIdentities] = getTrialFrameInds(minVel, obsCenter, obsPrePostBot, velPrePost, frameTimeStamps,...
    wheelPositions, wheelTimes, obsPositions, obsTimes, obsOnTimes, obsOffTimes);
frameIndsTop = getTrialFrameInds(minVel, obsCenter, obsPrePostTop, velPrePost, frameTimeStamps,...
    wheelPositions, wheelTimes, obsPositions, obsTimes, obsOnTimes, obsOffTimes);


vid = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4']);
showLocations(vid, frameIndsBot, potentialLocationsBot, ...
    locations.locationsCorrected, locations.trialIdentities, 1, .025, anchorPtsBot, hsv(4));










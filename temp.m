
% perform paw tracking for multiple sessions in parallel

sessionDirs = uigetdir2([getenv('OBSDATADIR') 'sessions\'], 'select folders to analyze');
steps = {'potBot', 'bot'};
minVel = .4;

for i = 1:length(sessionDirs)
    
    % get locations
    nameInd = find(sessionDirs{i}=='\',1,'last');
    try
        getLocations(sessionDirs{i}(nameInd+1:end), steps, minVel, false);
    catch
        fprintf('FAILED TO ANALYZE %s\n', sessionDirs{i}(nameInd+1:end));
    end
    
end

%% test getLocations on single session

session = '180122_002';
steps = {'potBot', 'bot'};
minVel = .6;
getLocations(session, steps, minVel, false);

% show tracking for session
load([getenv('OBSDATADIR') 'sessions\' session '\tracking\potentialLocationsBot.mat'])
load([getenv('OBSDATADIR') 'sessions\' session '\tracking\locationsBot.mat'])
vidBot = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4']);
anchorPtsBot = {[0 0], [1 0], [1 1], [0 1]}; % LH, LF, RF, RH // each entry is x,y pair measured from top left corner
showLocations(vidBot, find(locations.isAnalyzed), ...
        potentialLocationsBot, locations.locationsRaw, 1, .02, anchorPtsBot, hsv(4));

%% correct tracking

% settings
session = '180122_001';
view = 'Bot';
frameDelay = .03;

outputFile = [getenv('OBSDATADIR') 'sessions\' session '\tracking\locations' view 'Corrected.mat'];
load([getenv('OBSDATADIR') 'sessions\' session '\tracking\locations' view '.mat'], 'locations')
vid = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\run' view '.mp4']);
frameInds = find(locations.isAnalyzed);
anchorPts = {[0 0], [1 0], [1 1], [0 1]}; % LH, LF, RF, RH

correctTracking(outputFile, vid, locations, frameInds, frameDelay, anchorPts);







% perform paw tracking for multiple sessions in parallel

sessionDirs = uigetdir2([getenv('OBSDATADIR') 'sessions\'], 'select folders to analyze');
views = {'bot'};
minVel = .6;
%%

for i = 1:length(sessionDirs)
    
    % get locations
    nameInd = find(sessionDirs{i}=='\',1,'last');
    try
        getLocations(sessionDirs{i}(nameInd+1:end), views, minVel, false);
    catch
        fprintf('FAILED TO ANALYZE %s\n', sessionDirs{i}(nameInd+1:end));
    end
    
end


%% test correctTracking

% settings
session = '180122_000';
view = 'Bot';

outputFile = [getenv('OBSDATADIR') 'sessions\' session '\tracking\locations' view '.mat'];
load([getenv('OBSDATADIR') 'sessions\' session '\tracking\locations' view '.mat'], 'locations')
vid = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\run' view '.mp4']);
frameInds = find(locations.isAnalyzed);
anchorPts = {[0 0], [1 0], [1 1], [0 1]}; % LH, LF, RF, RH

correctTracking(outputFile, vid, locations, frameInds, .04, anchorPts);


%% rando tests

session = '180122_000';
vid = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4']);
load([getenv('OBSDATADIR') 'sessions\' session '\tracking\potentialLocationsBot.mat'], 'potentialLocationsBot')
load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], 'frameTimeStamps')
anchorPts = {[0 0], [1 0], [1 1], [0 1]}; % LH, LF, RF, RH
colors = hsv(4);

locationsBot = getLocationsBot(potentialLocationsBot, anchorPts, frameTimeStamps, vid.Width, vid.Height);
save([getenv('OBSDATADIR') 'sessions\' session '\tracking\locationsBot.mat'], 'locationsBot');
showLocations(vid, find([locationsBot.isAnalyzed]), potentialLocationsBot, locationsBot.locationsRaw, 1, .02, anchorPts, colors);







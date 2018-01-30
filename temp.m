
% perform paw tracking for multiple sessions in parallel

sessionDirs = uigetdir2([getenv('OBSDATADIR') 'sessions\'], 'select folders to analyze');
views = {'bot'};
minVel = .4;
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

session = '180122_000';
view = 'Bot';
load([getenv('OBSDATADIR') 'sessions\' session '\tracking\locations' view '.mat'], ['locations' view])
vid = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\run' view '.mp4']);
frameInds = find(locationsBot.isAnalyzed);
anchorPts = {[0 0], [1 0], [1 1], [0 1]}; % LH, LF, RF, RH

correctTracking(vid, eval(['locations' view]), frameInds, .04, anchorPts);








% function preparePoseRegressionData

% temp
sessions = {'180122_000', '180122_001', '180122_002'};
totalEgs = 500;

% settings
writeDir = 'C:\Users\rick\Desktop\trainingExamples\posRegression\';



% concatinate all labeled data sets
sessionInds = []; % stores the session identity for each saved location
locationsAll = [];

for i = 1:length(sessions)
    
    % get labeled locations for single session
    load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\tracking\locationsBotCorrected.mat'], 'locations');
    
    % remove nan entries
    locations = locations.locationsCorrected(locations.isAnalyzed,:,:);
    
    % store
    locationsAll = cat(1, locationsAll, locations);
    sessionInds = [sessionInds i*ones(1,size(locationsAll,1))];
    
end

locations = locationsAll; clear locationsAll;

% select random frames
locationInds = randperm(size(locations,1), totalEgs);
locationInds = sort(locationInds);











function getLocationsBot(showTracking)


% load tracking data
load('C:\Users\rick\Google Drive\columbia\obstacleData\svm\trackedData\potentialLocationsBot.mat', 'potentialLocationsBot')

% user settings
vidFile = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\testVideo\runBot.mp4';
dataDir = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\trackedData\';
objectNum = 4;
anchorPts = {[0 0], [0 1], [1 0], [1 1]}; % RH, LH, RF, LF (x, y)
maxDistanceX = .55;
maxVel = 35 / .004; % pixels / sec
minScore = 1.5; % location scores lower than minScores are set to zero (this way an object preferes occlusion to being assigned to a crummy location)

unaryWeight = 1.5;
pairwiseWeight = 1;
scoreWeight = 0;


% initializations
labels = nan(length(potentialLocationsBot), objectNum);
vid = VideoReader(vidFile);
frameTimes = 0:.004:.004*(vid.NumberOfFrames-1); % temp, these are fake timestamps! oh shit!



% set starting frame locations
unaries = nan(objectNum, length(potentialLocationsBot(1).x));

for i = 1:objectNum
    
    unaries(i,:) = getUnaryPotentials(potentialLocationsBot(1).x, potentialLocationsBot(1).y,...
                                      vid.Width, vid.Height,...
                                      anchorPts{i}(1), anchorPts{i}(2), maxDistanceX);
end


labels(1,:) = getBestLabels(unaries, objectNum, [0 0 0 0]);
locationsBot = struct();
locationsBot.x = nan(vid.NumberOfFrames, objectNum);
locationsBot.z = nan(vid.NumberOfFrames, objectNum);

locationsBot.x(1,:) = potentialLocationsBot(1).x(labels(1,:));
locationsBot.z(1,:) = potentialLocationsBot(1).y(labels(1,:));


% iterate through remaining frames

for i = 2:vid.NumberOfFrames
    
    % get unary and pairwise potentials
    unaries = nan(objectNum, length(potentialLocationsBot(i).x));
    pairwise = nan(objectNum, length(potentialLocationsBot(i).x));
    wasOccluded = zeros(1, objectNum);
    
    for j = 1:objectNum
        
        % unary
        unaries(j,:) = getUnaryPotentials(potentialLocationsBot(i).x, potentialLocationsBot(i).y,...
            vid.Width, vid.Height, anchorPts{j}(1), anchorPts{j}(2), maxDistanceX);
        
        % pairwise
        % get ind of last detection frame for object j
        if ~isnan(labels(i-1, j)) 
            prevFrame = i-1;
        else
            prevFrame = find(~isnan(labels(:,j)), 1, 'last');
            wasOccluded(j) = 1;
        end
        
        % get label at last dection frame
        prevLabel = labels(prevFrame, j);
        
        pairwise(j,:) = getPairwisePotentials(potentialLocationsBot(i).x, potentialLocationsBot(i).y,...
            potentialLocationsBot(prevFrame).x(prevLabel), potentialLocationsBot(prevFrame).y(prevLabel),...
            frameTimes(i)-frameTimes(prevFrame), maxVel);
        
    end
    
    % get track scores
    trackScores = repmat(potentialLocationsBot(i).scores', objectNum, 1);
    
    % find best labels
    scores = unaryWeight.*unaries + pairwiseWeight.*pairwise + scoreWeight.*trackScores;
    scores(unaries==0 | pairwise==0) = 0;
    scores(scores<minScore) = 0;
    labels(i,:) = getBestLabels(scores, objectNum, wasOccluded);
    
    % only keep labeled locations
    for j = 1:objectNum
        if isnan(labels(i,j))
            locationsBot.x(i,j) = nan;
            locationsBot.z(i,j) = nan;
        else
            locationsBot.x(i,j) = potentialLocationsBot(i).x(labels(i,j));
            locationsBot.z(i,j) = potentialLocationsBot(i).y(labels(i,j));
        end
    end
end

save([dataDir 'locationsBot.mat'], 'locationsBot');

% show tracking
if showTracking
    startFrame = 1;
    showPotentialLocations = true;
    showLocations(vid, potentialLocationsBot, locationsBot, showPotentialLocations, .04, anchorPts, startFrame);
end










% load tracking data
load('C:\Users\rick\Google Drive\columbia\obstacleData\svm\trackedData\trackedBotAll.mat', 'locationsBotAll')

% user settings
vidFile = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\testVideo\runBot.mp4';
objectNum = 4;
anchorPts = {[0 0], [0 1], [1 0], [1 1]}; % RH, LH, RF, LF (x, y)
maxDistanceX = .6;
maxVel = 35 / .004; % pixels / sec

unaryWeight = 1.5;
pairwiseWeight = 1;
scoreWeight = 0;


% initializations
labels = nan(length(locationsBotAll), objectNum);
vid = VideoReader(vidFile);
frameTimes = 0:.004:.004*(vid.NumberOfFrames-1); % temp, these are fake timestamps! oh shit!



% set starting frame locations
unaries = nan(objectNum, length(locationsBotAll(1).x));

for i = 1:objectNum
    
    unaries(i,:) = getUnaryPotentials(locationsBotAll(1).x, locationsBotAll(1).y,...
                                      vid.Width, vid.Height,...
                                      anchorPts{i}(1), anchorPts{i}(2), maxDistanceX);
end

labels(1,:) = getBestLabels(unaries, objectNum);
locationsBot = struct();
locationsBot(1).x = locationsBotAll(1).x(labels(1,:));
locationsBot(1).y = locationsBotAll(1).y(labels(1,:));


% iterate through remaining frames

for i = 2:vid.NumberOfFrames
    
    % get unary and pairwise potentials
    unaries = nan(objectNum, length(locationsBotAll(i).x));
    pairwise = ones(objectNum, length(locationsBotAll(i).x));
    isOccluded = zeros(1, objectNum);
    
    for j = 1:objectNum
        
        % unary
        unaries(j,:) = getUnaryPotentials(locationsBotAll(i).x, locationsBotAll(i).y,...
            vid.Width, vid.Height, anchorPts{j}(1), anchorPts{j}(2), maxDistanceX);
        
        % pairwise
        % get ind of last detection frame for object j
        if ~isnan(labels(i-1, j)) 
            prevFrame = i-1;
        else
            prevFrame = find(~isnan(labels(:,j)), 1, 'last');
            isOccluded(j) = 1;
        end
        
        % get label at last dection frame
        prevLabel = labels(prevFrame, j);
        
        pairwise(j,:) = getPairwisePotentials(locationsBotAll(i).x, locationsBotAll(i).y,...
            locationsBotAll(prevFrame).x(prevLabel), locationsBotAll(prevFrame).y(prevLabel),...
            frameTimes(i)-frameTimes(prevFrame), maxVel);
        
    end
    
    % get track scores
    trackScores = repmat(locationsBotAll(i).scores', objectNum, 1);
    
    % find best labels
    scores = unaryWeight.*unaries + pairwiseWeight.*pairwise + scoreWeight.*trackScores;
    scores(unaries==0 | pairwise==0) = 0;
    labels(i,:) = getBestLabels(scores, objectNum, isOccluded);
    
    % only keep labelled locations
    for j = 1:objectNum
        if isnan(labels(i,j))
            locationsBot(i).x(j) = nan;
            locationsBot(i).y(j) = nan;
        else
            locationsBot(i).x(j) = locationsBotAll(i).x(labels(i,j));
            locationsBot(i).y(j) = locationsBotAll(i).y(labels(i,j));
        end
    end
end

save([dataDir 'trackedBot.mat'], 'locationsBot');

%% show tracking
startFrame = 1;
showTracking(vid, locationsBotAll, labels, .04, anchorPts, startFrame);






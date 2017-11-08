

% load tracking data
load('C:\Users\rick\Google Drive\columbia\obstacleData\svm\trackedData\trackedBot.mat', 'locations')

% user settings
vidFile = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\testVideo\runBot.mp4';
objectNum = 4;
anchorPts = {[0 0], [0 1], [1 0], [1 1]}; % RH, LH, RF, LF (x, y)
maxDistanceX = .8;

unaryWeight = 1;
pairwiseWeight = 1;
scoreWeight = 1;


% initializations
labels = nan(length(locations), objectNum);
vid = VideoReader(vidFile);
frameTimes = 0:.004:.004*(vid.NumberOfFrames-1); % temp, these are fake timestamps! oh shit!



% set starting frame locations
unaries = nan(objectNum, length(locations(1).x));

for i = 1:objectNum
    
    unaries(i,:) = getUnaryPotentials(locations(1).x, locations(1).y,...
                                      vid.Width, vid.Height,...
                                      anchorPts{i}(1), anchorPts{i}(2), maxDistanceX);
end

labels(1,:) = getBestLabels(unaries, objectNum);



%% iterate through remaining frames

for i = 2:vid.NumberOfFrames
    
    % get unary and pairwise potentials
    unaries = nan(objectNum, length(locations(i).x));
    pairwise = ones(objectNum, length(locations(i).x));
    
    for j = 1:objectNum
        
        % unary
        unaries(j,:) = getUnaryPotentials(locations(i).x, locations(i).y,...
            vid.Width, vid.Height, anchorPts{j}(1), anchorPts{j}(2), maxDistanceX);
        
        % pairwise
        prevInd = 1; % !!! get last detected ind for object
        unaries(j,:) = getPairwisePotentials(locations(i).x, locations(i).y,...
            locations(prevInd).x, locations(prevInd).y,...
            frameTimes(i)-frameTimes(prevInd), maxVel);
        
    end
    
    
    for j = 1:objectNum
    
    % get track scores
    trackScores = repmat(locations(i).scores', objectNum, 1);
    
    % find best labels
    scores = unaryWeight.*unaries + pairwiseWeight.*pairwise + scoreWeight.*trackScores;
    labels(i,:) = getBestLabels(scores, objectNum);
end












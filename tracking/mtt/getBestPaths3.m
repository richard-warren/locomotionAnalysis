% to try: make frame skipping inferface // smarter reappearance function?

session = 'markerTest2';

% load tracking data
load([getenv('OBSDATADIR') 'sessions/' session '/tracking/potentialLocationsBot.mat'], 'potentialLocationsBot');

% user settings
unaryWeight = 1;
pairwiseWeight = 0;
occludedWeight = .001;
occlusionGridSpacing = 50;
maxDistanceX = .65;
maxDistanceY = .65;
maxVelocity = 30;
anchorPts = {[0 0], [0 1], [1 0], [1 1]}; % LH, RH, LF, RF (x, y)
paws = [2];

% initializations
frameInds = 34757:35085; % temp
vid = VideoReader([getenv('OBSDATADIR') 'sessions/' session '/runBot.mp4']);
[gridX, gridY] = meshgrid(1:occlusionGridSpacing:vid.Width,...
                          1:occlusionGridSpacing:vid.Height);
gridPts = [gridX(:), gridY(:)];
numOccluded = size(gridPts,1);
locations = potentialLocationsBot;
locationsBot.x = nan(length(potentialLocationsBot), 4);
locationsBot.y = nan(length(potentialLocationsBot), 4);


%% viterbi

tic
locationScores = cell(length(frameInds), 1);
locationTraceBacks = cell(length(frameInds), 1);

for j = paws
    
    % forward
    disp('forward!')
    for i = 1:length(frameInds)
    
        currentFrame = frameInds(i);
        currentNum = length(locations(currentFrame).x);
        
        % get unary potentials
        unaries = getUnaryPotentials(locations(currentFrame).x, locations(currentFrame).y, vid.Width, vid.Height,...
            anchorPts{j}(1), anchorPts{j}(2), maxDistanceX, maxDistanceY);
        unaries = [unaries * unaryWeight; zeros(numOccluded,1)];
        
        % get pairwise potentials
        if i>1
            pairwise = getPairwisePotentials([locations(currentFrame).x, locations(currentFrame).y],...
                [locations(currentFrame-1).x, locations(currentFrame-1).y],...
                maxVelocity, pairwiseWeight, occludedWeight, gridPts);
        else
            pairwise = zeros(currentNum + numOccluded, 1);
        end
        
        % !!! need to make invalid transitions impossible
        
        % get best score for each potential location, as well as previous location that led to that score        
        scores = pairwise + repmat(unaries, 1, size(pairwise,2));
        [locationScores{i}, locationTraceBacks{i}] = max(scores, [], 2);
        % !!! need to multiply current scores by previous best scores...
    end
    
    
    % backward
    disp('backward!')
    [~, prevInd] = max(locationScores{end});
        
    for i = fliplr(1:length(frameInds)-1)
        
        if prevInd <= length(potentialLocationsBot(frameInds(i)).x)
            locationsBot.x(frameInds(i), j) = potentialLocationsBot(frameInds(i)).x(prevInd);
            locationsBot.y(frameInds(i), j) = potentialLocationsBot(frameInds(i)).y(prevInd);
        end
        
        prevInd = locationTraceBacks{i}(prevInd);
    end
end


% visualize tracking
showLocations(vid, frameInds, potentialLocationsBot, locationsBot, true, .02, anchorPts);





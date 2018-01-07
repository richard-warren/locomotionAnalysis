% to try: make frame skipping inferface // smarter reappearance function?

session = 'markerTest2';

% load tracking data
load([getenv('OBSDATADIR') 'sessions/' session '/tracking/potentialLocationsBot.mat'], 'potentialLocationsBot');

% user settings
unaryWeight = 1;
pairwiseWeight = .1;
occludedWeight = .001;
occlusionGridSpacing = 30;
maxDistanceX = .65;
maxDistanceY = .65;
maxVelocity = 30;
anchorPts = {[0 0], [0 1], [1 0], [1 1]}; % LH, RH, LF, RF (x, y)

% initializations
frameInds = 34757:35085; % temp
vid = VideoReader([getenv('OBSDATADIR') 'sessions/' session '/runBot.mp4']);
[gridX, gridY] = meshgrid(1:occlusionGridSpacing:vid.Width,...
                          1:occlusionGridSpacing:vid.Height);
gridPts = [gridX(:), gridY(:)];
numOccluded = size(gridPts,1);
potentialLocations = potentialLocationsBot;


%% viterbi

tic
locations.x = nan(length(potentialLocationsBot), 4);
locations.y = nan(length(potentialLocationsBot), 4);
locationScores = cell(length(frameInds), 1);
locationTraceBacks = cell(length(frameInds), 1);

paws = [2];
close all; figure; im = imagesc(randn(10));

for j = paws
    
    % forward
    disp('forward!')
    for i = 1:length(frameInds)
    
        currentFrame = frameInds(i);
        currentNum = length(potentialLocations(currentFrame).x);
        
        % get unary potentials
        [unaries, invalidPositions] = getUnaryPotentials(potentialLocations(currentFrame).x, potentialLocations(currentFrame).y, vid.Width, vid.Height,...
            anchorPts{j}(1), anchorPts{j}(2), maxDistanceX, maxDistanceY);
        unaries = [unaries * unaryWeight; zeros(numOccluded,1)];
        
        % get pairwise potentials
        if i>1
            [pairwise, invalidTransitions] = getPairwisePotentials([potentialLocations(currentFrame).x, potentialLocations(currentFrame).y],...
                [potentialLocations(currentFrame-1).x, potentialLocations(currentFrame-1).y],...
                maxVelocity, pairwiseWeight, occludedWeight, gridPts);
            prevNum = length(potentialLocations(currentFrame-1).x);
        else
            pairwise = ones(currentNum + numOccluded, 1);
            prevNum = 1;
            invalidTransitions = true(size(pairwise));
        end
        
        % get best score for each potential location, as well as previous location that led to that score        
        scores = pairwise;
        scores(:,1:prevNum) = scores(:,1:prevNum) + repmat(unaries, 1, prevNum);
        scores(:,1:prevNum) = scores(:,1:prevNum) + repmat(unaries, 1, prevNum);
        scores(invalidTransitions) = 0;
        scores(invalidPositions, 1:prevNum) = 0;
        
        
        % !!! adjust getPairewisePotentials to return invalid transitions
        
        [locationScores{i}, locationTraceBacks{i}] = max(scores, [], 2);
        if i>1 % multiply current scores by previous best scores
            locationScores{i} = locationScores{i} .* locationScores{i-1}(locationTraceBacks{i});
%             locationScores{i} = locationScores{i} / nansum(locationScores{i});
        end
        
%         set(im, 'CData', scores(1:min(10,size(scores,1)),1:min(10,size(scores,2))))
%         pause(1)
        
    end
    
    
    % backward
    disp('backward!')
    [~, currentInd] = max(locationScores{end});
    
    for i = fliplr(1:length(frameInds))
        
%         [~, maxInd] = max(locationScores{i});
%         if maxInd <= length(potentialLocationsBot(frameInds(i)).x) % if it is not occluded
%             locations.x(frameInds(i),j) = potentialLocations(frameInds(i)).x(maxInd);
%             locations.y(frameInds(i),j) = potentialLocations(frameInds(i)).y(maxInd);
%         end
        
        if currentInd <= length(potentialLocationsBot(frameInds(i)).x) % if it is not occluded
            locations.x(frameInds(i), j) = potentialLocations(frameInds(i)).x(currentInd);
            locations.y(frameInds(i), j) = potentialLocations(frameInds(i)).y(currentInd);
        end
%         
%         currentInd = locationTraceBacks{i}(currentInd);
    end
end


% visualize tracking
showLocations(vid, frameInds, potentialLocations, locations, true, .02, anchorPts);





% to try: make frame skipping inferface // smarter reappearance function?

session = 'markerTest2';


% user settings
frameInds = 36120:36499; % temp

unaryWeight = 1;
pairwiseWeight = .1;
occludedWeight = .1;
occlusionGridSpacing = 50;
maxDistanceX = .65;
maxDistanceY = .65;
maxVelocity = 30;
anchorPts = {[0 0], [0 1], [1 0], [1 1]}; % LH, RH, LF, RF (x, y)

% initializations
load([getenv('OBSDATADIR') 'sessions/' session '/tracking/potentialLocationsBot.mat'], 'potentialLocationsBot');
potentialLocations = potentialLocationsBot; clear potentialLocationsBot
vid = VideoReader([getenv('OBSDATADIR') 'sessions/' session '/runBot.mp4']);


[gridX, gridY] = meshgrid(1:occlusionGridSpacing:vid.Width,...
                          1:occlusionGridSpacing:vid.Height);
gridPts = [gridX(:), gridY(:)];
numOccluded = size(gridPts,1);


%% viterbi

locations.x = nan(length(potentialLocations), 4);
locations.y = nan(length(potentialLocations), 4);
locationScores = cell(length(frameInds), 1);
locationTraceBacks = cell(length(frameInds), 1);

paws = 1:2;
% close all; figure; im = imagesc(randn(10));

for j = paws
    
    % forward
    disp('forward!')
    
    for i = 1:length(frameInds)
    
        currentFrame = frameInds(i);
        currentNum = length(potentialLocations(currentFrame).x);
        
        % get unary potentials
        [unaries, invalidPositions] = getUnaryPotentials(potentialLocations(currentFrame).x, potentialLocations(currentFrame).y,...
            vid.Width, vid.Height, anchorPts{j}(1), anchorPts{j}(2), maxDistanceX, maxDistanceY);
        unaries = [unaries * unaryWeight; zeros(numOccluded,1)];
        
        % get pairwise potentials
        if i>1 
            [pairwise, invalidTransitions] = getPairwisePotentialsViterbi([potentialLocations(currentFrame).x, potentialLocations(currentFrame).y],...
                [potentialLocations(currentFrame-1).x, potentialLocations(currentFrame-1).y],...
                maxVelocity, pairwiseWeight, occludedWeight, gridPts);
            prevNum = length(potentialLocations(currentFrame-1).x);
        else
            pairwise = [ones(currentNum, 1); zeros(numOccluded, 1)];
            prevNum = 1;
            invalidTransitions = false(currentNum,1);
        end
        
        % compute scores
        scores = pairwise;
        scores(:,1:prevNum) = scores(:,1:prevNum) + repmat(unaries, 1, prevNum);
        scores(1:currentNum, 1:prevNum) = scores(1:currentNum, 1:prevNum) .* ~invalidTransitions;
        scores(1:currentNum, 1:prevNum) = scores(1:currentNum, 1:prevNum) .* ~repmat(invalidPositions,1,prevNum);
        
        [bestScores, locationTraceBacks{i}] = max(scores, [], 2);
        if i>1; bestScores = bestScores .* locationScores{i-1}(locationTraceBacks{i}); end
        locationScores{i} = bestScores;
        if max(bestScores)==0; fprintf('zero max: %i\n', j); end
                
%         set(im, 'CData', scores(1:min(10,size(scores,1)),1:min(10,size(scores,2))))
%         pause(1)
        
    end
    
    
    % backward
    disp('backward!')
    [~, currentInd] = max(locationScores{end});
    
    for i = fliplr(1:length(frameInds))
        
%         [~, maxInd] = max(locationScores{i});
%         if maxInd <= length(potentialLocations(frameInds(i)).x) % if it is not occluded
%             locations.x(frameInds(i),j) = potentialLocations(frameInds(i)).x(maxInd);
%             locations.y(frameInds(i),j) = potentialLocations(frameInds(i)).y(maxInd);
%         end
        
        if currentInd <= length(potentialLocations(frameInds(i)).x) % if it is not occluded
            locations.x(frameInds(i),j) = potentialLocations(frameInds(i)).x(currentInd);
            locations.y(frameInds(i),j) = potentialLocations(frameInds(i)).y(currentInd);
        end
        
        currentInd = locationTraceBacks{i}(currentInd);
    end
end


% visualize tracking
showLocations(vid, frameInds, potentialLocations, locations, true, .02, anchorPts);





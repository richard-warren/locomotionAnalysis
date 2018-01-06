% to try: make frame skipping inferface // smarter reappearance function?

session = 'markerTest2';

% load tracking data
load([getenv('OBSDATADIR') 'sessions/' session '/tracking/potentialLocationsBot.mat'], 'potentialLocationsBot');

% user settings
unaryWeight = 1;
pairwiseWeight = 1;
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

tic

%% viterbi forward

locationScores = cell(length(frameInds), 1);
locationTraceBacks = cell(length(frameInds), 1);

for i = 1:length(frameInds)
    
    currentFrame = frameInds(i);

    for j = paws
        
        % get unary potentials
        unaries = getUnaryPotentials(locations(currentFrame).x, locations(currentFrame).y, vid.Width, vid.Height,...
            anchorPts{j}(1), anchorPts{j}(2), maxDistanceX, maxDistanceY);
        unaires = unaries * unaryWeight;
        
        % get pairwise potentials
        if i>1
            pairwise = getPairwisePotentials([locations(currentFrame).x, locations(currentFrame).y],...
                [locations(currentFrame-1).x, locations(currentFrame-1).y],...
                maxVelocity, pairwiseWeight, occludedWeight, gridPts);
        else
            pairwise = zeros(length(locations(currentFrame).x));
        end
        
        % get best score for each potential location, as well as previous location that led to that score
        if i>1
            currentNum = length(unaries);
            prevNum = length(locationScores(getPairwisePotentials-1).x);
            
            scores = pairwise;
            scores(1:currentNum, 1:prevNum) = scores(1:currentNum, 1:prevNum) + ...
                                            + repmat(unaries, 1, prevNum));
            locationScores{i} = scores;
            locationTraceBacks{i};
        else
        end
    
    
    
    
    
    end
        
    
    
    % get pairwise potentials
    
end



%% find most probable paths!

% initializations
objectNum = size(unary{1}, 1);
labels = nan(length(unary), objectNum);
nodeScores = cell(1, length(unary));
backPointers = cell(1, length(unary)-1);
pathScores = nan(1,4);
nodeScores{1} = unary{1};


for i = 1:objectNum
    
    for j = 2:length(unary)
        
        % forward
        
        % initialize containers for first object only
        if i==1
            nodeScores{j}   = nan(objectNum, size(pairwise{j-1},1));
            backPointers{j-1} = nan(objectNum, size(pairwise{j-1},1));
        end
        
        previousScores = nodeScores{j-1}(i,:);
        currentUnary   = unary{j}(i,:)';
        
        allTransitionScores = repmat(previousScores, length(currentUnary), 1) + ...
                              repmat(currentUnary, 1, length(previousScores)) + ...
                              pairwise{j-1};
        
        % make invalid transitions impossible                      
        invalidInds = (pairwise{j-1}==0) |...
                      (repmat(previousScores, length(currentUnary), 1)==0);          
        allTransitionScores(invalidInds) = 0;
        
        % make transitions to locations with 0 unary potential impossible
        inds = find(currentUnary==0);
        allTransitionScores(inds<=length(locations(j).x), 1:length(locations(j-1).x)) = 0;
        
        % make transitions to previously occupied state impossible
        if i>1
            labelInds = labels(j-1,1:i-1);
            labelInds = labelInds( labelInds <= length(locations(j-1).x) ); % occluded locations can be multiply occupied
            allTransitionScores(:, labelInds) = 0;
        end
        
        
        [nodeScores{j}(i,:), backPointers{j-1}(i,:)] = max(allTransitionScores, [], 2);
%         backPointers{j-1}(i,:) = backPointers{j-1}(i,:) .* double(nodeScores{j}(i,:)>0);
    end
    
    % backward
    [pathScores(i), labels(end,i)] = max(nodeScores{end}(i,:));
    
    for j = fliplr(1:length(unary)-1) 
        labels(j,i) = backPointers{j}(i, labels(j+1,i));
    end
end


toc

%% visualize tracking
showTracking(VideoReader(vid), locations, labels, gridPts, vidDelay, paws);





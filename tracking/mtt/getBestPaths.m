% to try: make frame skipping inferface // smarter reappearance function?


% load tracking data
load('C:\Users\rick\Google Drive\columbia\obstacleData\svm\trackedData\tracked.mat', 'locations')

% user settings
vidFile = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\testVideo\botTest.mp4';

unaryWeight = 2;
pairwiseWeight = .1;
occludedWeight = .01;
occlusionGridSpacing = 20;

maxVelocity = 30;
paws = 1:4;
vidDelay = .02;

% initializations
frameHeight = 242; % !!! hacky temp
frameWidth = 398;
[gridX, gridY] = meshgrid(1:occlusionGridSpacing:frameWidth,...
                          1:occlusionGridSpacing:frameHeight);
gridPts = [gridX(:), gridY(:)];
numOccluded = size(gridPts,1);


% compute unary potentials
unary = cell(1,length(locations));

for i = 1:length(locations)

    % ensure first unaries are not zero (this will cause all paths to be zero for all time!)
    if (i==1 && unaryWeight==0)
        unaryWeightTemp = 1;
    else
        unaryWeightTemp = unaryWeight;
    end

    frameUnaries = nan(4, length(locations(i).x) + numOccluded);
    
    frameUnaries(1,1:end) = getUnaryPotentials(locations(i).x, locations(i).y, frameWidth, frameHeight, 0, .2, numOccluded, unaryWeightTemp); % RH
    frameUnaries(2,1:end) = getUnaryPotentials(locations(i).x, locations(i).y, frameWidth, frameHeight, 0, 1, numOccluded, unaryWeightTemp); % RF
    frameUnaries(3,1:end) = getUnaryPotentials(locations(i).x, locations(i).y, frameWidth, frameHeight, 1, .2, numOccluded, unaryWeightTemp); % LH
    frameUnaries(4,1:end) = getUnaryPotentials(locations(i).x, locations(i).y, frameWidth, frameHeight, 1, 1, numOccluded, unaryWeightTemp); % LF
    
    
    unary{i} = frameUnaries;
    
end



% compute pairwise potentials
pairwise = cell(1,length(locations)-1);

for i = 1:length(locations)-1
    
    pairwise{i} = getPairwisePotentials([locations(i+1).x, locations(i+1).y], [locations(i).x, locations(i).y],...
                                        maxVelocity, pairwiseWeight, occludedWeight, gridPts);
    
end


%% attempt at match2nd

unaryFlipped = cellfun(@(x) x', unary, 'un', 0);
pairwiseFlipped = cellfun(@(x) x', pairwise, 'un', 0);

labels = match2nd(unaryFlipped, pairwise, [], numOccluded, 0);
showTracking(VideoReader(vidFile), locations, labels, gridPts, vidDelay, paws);


%% find most probable paths!

objectNum = size(unary{1}, 1);
labels = nan(length(unary), objectNum);
nodeScores = cell(1, length(unary));
backPointers = cell(1, length(unary)-1);

% initializations

nodeScores{1} = unary{1};


for i = 2:length(unary)
    
    nodeScores{i}   = nan(objectNum, size(pairwise{i-1},1));
    backPointers{i-1} = nan(objectNum, size(pairwise{i-1},1));
    
    for j = 1:objectNum
        
        previousScores = nodeScores{i-1}(j,:);
        currentUnary   = unary{i}(j,:)';
        
        allTransitionScores = repmat(previousScores, length(currentUnary), 1) + ...
                              repmat(currentUnary, 1, length(previousScores)) + ...
                              pairwise{i-1};
                              
        invalidInds = (pairwise{i-1}==0) | (repmat(previousScores, length(currentUnary), 1)==0);
        allTransitionScores(invalidInds) = 0; % make invalid transitions impossible
        
        [nodeScores{i}(j,:), backPointers{i-1}(j,:)] = max(allTransitionScores, [], 2);
        backPointers{i-1}(j,:) = backPointers{i-1}(j,:);% .* double(nodeScores{i}(j,:)>0);
    end
end


% back trace

[pathScores, labels(end,:)] = max(nodeScores{end}, [], 2);
%
for i = fliplr(1:length(unary)-1)
    for j = 1:objectNum
        
        labels(i,j) = backPointers{i}(j, labels(i+1,j));

    end
end




% VISUALIZE TRACKING
showTracking(VideoReader(vidFile), locations, labels, gridPts, vidDelay, paws);





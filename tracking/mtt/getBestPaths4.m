% to try: make frame skipping inferface // smarter reappearance function?

session = 'markerTest2';


% user settings
frameInds = 36120:36499; % temp

unaryWeight = 1;
pairwiseWeight = .1;
occludedWeight = .01;
occlusionGridSpacing = 30;
maxDistanceX = .65;
maxDistanceY = .65;
maxVelocity = 50;
anchorPts = {[0 0], [0 1], [1 0], [1 1]}; % LH, RH, LF, RF (x, y)

% initializations
load([getenv('OBSDATADIR') 'sessions/' session '/tracking/potentialLocationsBot.mat'], 'potentialLocationsBot');
potentialLocations = potentialLocationsBot; clear potentialLocationsBot
vid = VideoReader([getenv('OBSDATADIR') 'sessions/' session '/runBot.mp4']);


[gridX, gridY] = meshgrid(1:occlusionGridSpacing:vid.Width,...
                          1:occlusionGridSpacing:vid.Height);
gridPts = [gridX(:), gridY(:)];
numOccluded = size(gridPts,1);


% viterbi

locations.x = nan(length(potentialLocations), 4);
locations.y = nan(length(potentialLocations), 4);
locationScores = cell(length(frameInds), 1);
locationTraceBacks = cell(length(frameInds), 1);

paws = 1;
paws = paws(randperm(length(paws)));
% close all; figure; im = imagesc(randn(10));
unaries = cell(1, length(frameInds));
pairwise = cell(1, length(frameInds)-1);


    
% get unary and pairwise potentials

for i = 1:length(frameInds)

    currentFrame = frameInds(i);
    currentNum = length(potentialLocations(currentFrame).x);
    unaries{i} = zeros(currentNum+numOccluded, length(paws));

    % get unary potentials
    for j = paws
        unaries{i}(1:currentNum, j) = getUnaryPotentials(potentialLocations(currentFrame).x, potentialLocations(currentFrame).y,...
            vid.Width, vid.Height, anchorPts{j}(1), anchorPts{j}(2), maxDistanceX, maxDistanceY);
        unaries{i}(1:currentNum, j) = unaries{i}(1:currentNum, j) * unaryWeight;
    end

    % get pairwise potentials
    if i>1
        pairwise{i-1} = getPairwisePotentialsViterbi([potentialLocations(currentFrame).x, potentialLocations(currentFrame).y],...
            [potentialLocations(currentFrame-1).x, potentialLocations(currentFrame-1).y],...
            maxVelocity, pairwiseWeight, occludedWeight, gridPts);
        pairwise{i-1} = sparse(pairwise{i-1});
    end
end



% russell et al attempt

labels = match2nd(unaries, pairwise, [], numOccluded, 0);


for i = 1:length(frameInds)
    for j = paws
        if labels(j,i) <= length(potentialLocations(frameInds(i)).x) % if it is not occluded
            locations.x(frameInds(i),j) = potentialLocations(frameInds(i)).x(labels(j,i));
            locations.y(frameInds(i),j) = potentialLocations(frameInds(i)).y(labels(j,i));
        end
    end
end

showLocations(vid, frameInds, potentialLocations, fixTracking(locations), true, .02, anchorPts);








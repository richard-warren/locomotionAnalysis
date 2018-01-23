
session = '180109_002';


% initializations
load([getenv('OBSDATADIR') 'sessions/' session '/tracking/potentialLocationsBot.mat'], 'potentialLocationsBot');
potentialLocations = potentialLocationsBot; clear potentialLocationsBot
vid = VideoReader([getenv('OBSDATADIR') 'sessions/' session '/runBot.mp4']);

%% settings
trial = 100;
paws = 1:4;
unaryWeight = 1;
pairwiseWeight = .1;
occludedWeight = .001;
occlusionGridSpacing = 75;
maxDistanceX = .5;
maxDistanceY = .65;
maxVelocity = 50;
anchorPts = {[0 0], [0 1], [1 0], [1 1]}; % LH, RH, LF, RF (x, y)



% initializations
frameInds = find(~isnan(obsPixPositions));
trialStartInds = find(diff(frameInds)>250);
frameInds = frameInds(trialStartInds(trial)+1 : trialStartInds(trial+1));

[gridX, gridY] = meshgrid(1:occlusionGridSpacing:vid.Width,...
                          1:occlusionGridSpacing:vid.Height);
gridPts = [gridX(:), gridY(:)];
numOccluded = size(gridPts,1);

locations.x = nan(length(potentialLocations), 4);
locations.y = nan(length(potentialLocations), 4);
unaries = cell(1, length(frameInds));
pairwise = cell(1, length(frameInds)-1);



% get unary and pairwise potentials

for i = 1:length(frameInds)

    currentFrame = frameInds(i);
    currentNum = length(potentialLocations(currentFrame).x);
    unaries{i} = zeros(currentNum+numOccluded, length(paws));

    % get unary potentials
    for j = paws
        unaries{i}(1:currentNum, j) = getUnaryPotentialsViterbi(potentialLocations(currentFrame).x, potentialLocations(currentFrame).y,...
            vid.Width, vid.Height, anchorPts{j}(1), anchorPts{j}(2), maxDistanceX, maxDistanceY);
        unaries{i}(1:currentNum, j) = unaries{i}(1:currentNum, j) * unaryWeight;
    end
    
%     if frameInds(i)==36325; keyboard; end

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




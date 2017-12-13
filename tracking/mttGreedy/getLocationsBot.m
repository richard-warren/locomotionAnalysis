function locationsBot = getLocationsBot(potentialLocationsBot, frameTimeStamps, frameWid, frameHgt, frameInds)

% !!! need to document


% settings
objectNum = 4;
anchorPts = {[0 0], [0 1], [1 0], [1 1]}; % RH, LH, RF, LF (x, y)
maxDistanceX = .6;    % x coordinates can only be this far away from the x anchor points (expressed as percentage of frame width)
maxDistanceY = .55;
maxVel = 30 / .004;    % pixels / sec


unaryWeight = 1;
pairwiseWeight = 1;
minScore = .5 * (unaryWeight + pairwiseWeight); % location scores lower than minScores are set to zero (this way an object preferes occlusion to being assigned to a crummy location)
scoreWeight = 0;


% initializations
labels = nan(length(potentialLocationsBot), objectNum);


% set starting frame locations
unaries = nan(objectNum, length(potentialLocationsBot(frameInds(1)).x));

for i = 1:objectNum
    
    unaries(i,:) = getUnaryPotentials(potentialLocationsBot(frameInds(1)).x, potentialLocationsBot(frameInds(1)).y,...
                                      frameWid, frameHgt,...
                                      anchorPts{i}(1), anchorPts{i}(2), maxDistanceX, maxDistanceY);
end

labels(frameInds(1),:) = getBestLabels(unaries, objectNum, [0 0 0 0]);
locationsBot = struct();
locationsBot.x = nan(length(potentialLocationsBot), objectNum);
locationsBot.y = nan(length(potentialLocationsBot), objectNum);

locationsBot.x(frameInds(1),:) = potentialLocationsBot(frameInds(1)).x(labels(frameInds(1),:));
locationsBot.y(frameInds(1),:) = potentialLocationsBot(frameInds(1)).y(labels(frameInds(1),:));



% iterate through remaining frames

for i = frameInds(2:end)%length(potentialLocationsBot)
    
    % get unary and pairwise potentials
    unaries = nan(objectNum, length(potentialLocationsBot(i).x));
    pairwise = nan(objectNum, length(potentialLocationsBot(i).x));
    wasOccluded = zeros(1, objectNum);
    
    for j = 1:objectNum
        
        % unary
        unaries(j,:) = getUnaryPotentials(potentialLocationsBot(i).x, potentialLocationsBot(i).y,...
            frameWid, frameHgt, anchorPts{j}(1), anchorPts{j}(2), maxDistanceX, maxDistanceY);
        
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
            frameTimeStamps(i)-frameTimeStamps(prevFrame), maxVel);
        
    end
    
    % get track scores
    trackScores = repmat(potentialLocationsBot(i).scores', objectNum, 1);
    
    % find best labels
    scores = unaryWeight.*unaries + pairwiseWeight.*pairwise + scoreWeight.*trackScores;
    scores(unaries==0 | pairwise==0) = 0;
    scores(scores<minScore) = 0;
    
    % !!! this if/then is temporary // it would be better if getBestLabels handled empty values itself // NEED TO LOOK INTO WHAT HAPPENS WHEN THERE ARE FEWER POTENTIAL OBJECTS THAN LOCATIONS... ARE NANS RETURNED?
    if isempty(scores)
        labels(i,:) = nan;
    else
        labels(i,:) = getBestLabels(scores, objectNum, wasOccluded);
    end
    
    % only keep labeled locations
    for j = 1:objectNum
        if isnan(labels(i,j))
            locationsBot.x(i,j) = nan;
            locationsBot.y(i,j) = nan;
        else
            locationsBot.x(i,j) = potentialLocationsBot(i).x(labels(i,j));
            locationsBot.y(i,j) = potentialLocationsBot(i).y(labels(i,j));
        end
    end
    
    disp(i/length(potentialLocationsBot))
end






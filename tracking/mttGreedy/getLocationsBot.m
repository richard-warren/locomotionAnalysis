function locationsBot = getLocationsBot(potentialLocationsBot, frameTimeStamps, frameWid, frameHgt, frameInds)

% !!! need to document


% settings
objectNum = 4;
anchorPts = {[0 0], [0 1], [1 0], [1 1]}; % LH, RH, LF, RF (x, y)
maxDistanceX = .65;    % x coordinates can only be this far away from the x anchor points (expressed as percentage of frame width)
maxDistanceY = .65;
maxVel = 30 / .004;    % pixels / sec
pawCheckXProximity = 40; % if x locations of paws are within this value of one another, ensure forepaw is more medial than hind paw
pawCheckYDistance = 10;

unaryWeight = 2;
pairwiseWeight = 1;
minScore = .5 * (unaryWeight + pairwiseWeight); % location scores lower than minScores are set to zero (this way an object prefers occlusion to being assigned to a crummy location)
scoreWeight = 0;


% initializations
labels = nan(length(potentialLocationsBot), objectNum);
locationsBot.x = nan(length(potentialLocationsBot), objectNum);
locationsBot.y = nan(length(potentialLocationsBot), objectNum);
yMid = frameHgt / 2;
pawPairs = {[1,3], [2,4]}; % left paws (hind, fore), right paws (hind, fore)
pawSwaps = 0;


% iterate through remaining frames

for i = frameInds(1:end)
    
    % report progress
    disp(i/length(potentialLocationsBot))
    
    
    % get unary and pairwise potentials
    unaries = nan(objectNum, length(potentialLocationsBot(i).x));
    pairwise = ones(objectNum, length(potentialLocationsBot(i).x));
    wasOccluded = zeros(1, objectNum);
    
    for j = 1:objectNum
        
        % unary
        unaries(j,:) = getUnaryPotentials(potentialLocationsBot(i).x, potentialLocationsBot(i).y,...
            frameWid, frameHgt, anchorPts{j}(1), anchorPts{j}(2), maxDistanceX, maxDistanceY);
        
        % pairwise
        if i>frameInds(1)
            
            % get ind of last detection frame for object j
            if ~isnan(labels(i-1, j)) 
                prevFrame = i-1;
            else
                prevFrame = find(~isnan(labels(:,j)), 1, 'last');
                wasOccluded(j) = 1;
            end
            
            if ~isempty(prevFrame)
                
                % get label at last dection frame
                prevLabel = labels(prevFrame, j);

                % get pairwise scores
                pairwise(j,:) = getPairwisePotentials(potentialLocationsBot(i).x, potentialLocationsBot(i).y,...
                    potentialLocationsBot(prevFrame).x(prevLabel), potentialLocationsBot(prevFrame).y(prevLabel),...
                    frameTimeStamps(i)-frameTimeStamps(prevFrame), maxVel);
            end
        end
    end
    
    
    % get svm scores for potential locations
    svmScores = repmat(potentialLocationsBot(i).scores', objectNum, 1);
    
    
    % get final scores for each potential location for each paw
    scores = unaryWeight.*unaries + pairwiseWeight.*pairwise;% + scoreWeight.*svmScores;
    scores(unaries==0 | pairwise==0) = 0;
    scores(scores<minScore) = 0;

    
    % find best labels (find labels for all paws that maximizes overall score)
    if ~isempty(scores)
        labels(i,:) = getBestLabels(scores, objectNum, wasOccluded);
    end
    
    
    % if paws on same side (left or right) or both in the middle of frame, make sure the forepaw has the more medial label!
    for j = 1:length(pawPairs) % left paws (hind, fore), right paws (hind, fore)
        
        if ~any(isnan(labels(i, pawPairs{j}))) % if both paws on the side have labels
        
            xHind = potentialLocationsBot(i).x(labels(i, pawPairs{j}(1)));
            xFore = potentialLocationsBot(i).x(labels(i, pawPairs{j}(2)));

%             if xHind>overlapZone(1) && xHind<overlapZone(2) && ...
%                xFore>overlapZone(1) && xFore<overlapZone(2)
            if abs(xHind - xFore) < pawCheckXProximity % if paws are very close to eachother

                yHind = potentialLocationsBot(i).y(labels(i, pawPairs{j}(1)));
                yFore = potentialLocationsBot(i).y(labels(i, pawPairs{j}(2)));

                if (abs(yFore-yMid) - abs(yHind-yMid)) > pawCheckYDistance % if the y locations are separated by at least pawCheckYDistance pixels
                    labels(i, pawPairs{j}) = labels(i, fliplr(pawPairs{j}));
                    pawSwaps = pawSwaps + 1;
                    fprintf('paw swap #%i\n', pawSwaps)
                end

            end
        end
    end
    
    
    % only keep labeled locations
    for j = 1:objectNum
        if ~isnan(labels(i,j))
            locationsBot.x(i,j) = potentialLocationsBot(i).x(labels(i,j));
            locationsBot.y(i,j) = potentialLocationsBot(i).y(labels(i,j));
        end
    end
end

fprintf('total paw swaps: %i\n', pawSwaps)






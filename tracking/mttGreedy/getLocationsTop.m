function getLocationsTop(showTracking)

% !!! need to document and make not shitty

% load tracking data
load('C:\Users\rick\Google Drive\columbia\obstacleData\svm\trackedData\locationsBot.mat', 'locationsBot')
load('C:\Users\rick\Google Drive\columbia\obstacleData\svm\trackedData\potentialLocationsTop.mat', 'potentialLocationsTop')

% user settings
vidFile = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\testVideo\runTop.mp4';
dataDir = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\trackedData\';
objectNum = 4;
maxVel = 25 / .004; % pixels / sec
maxXDistance = 20;
xOccludeBuffer = 40;
minY = 70;

unariesWeight = 1;
pairwiseWeight = 1;
lownessWeight = 10;
scoreWeight = 0;


% initializations
labels = nan(length(potentialLocationsTop), objectNum);
locationsTop = struct();
vid = VideoReader(vidFile);
frameTimes = 0:.004:.004*(vid.NumberOfFrames-1); % temp, these are fake timestamps! oh shit!



% fill short periods of missing values in bottom paw tracking
locationsBot = fixTracking(locationsBot);



% iterate through remaining frames

for i = 1:vid.NumberOfFrames
    
    unaries = zeros(objectNum, length(potentialLocationsTop(i).x));
    pairwise = zeros(objectNum, length(potentialLocationsTop(i).x));
    valid = ones(objectNum, length(potentialLocationsTop(i).x));
    wasOccluded = zeros(1, objectNum); % keeps track of whether the object was occluded in the previous frame (used in getBestLabels)
    
    
    for j = 1:objectNum
        
        % check if paw is occluded
        occludedByBins = abs(locationsBot.x(i,j) - locationsBot.x(i,:)) < xOccludeBuffer & ...
                      (locationsBot.y(i,:) > locationsBot.y(i,j));
        if any(occludedByBins)
            valid(j,:) = 0;
        
        % if not occluded, compute scores for all potential locations
        else
            
            % get unary potentials
            xDistances = abs(locationsBot.x(i,j) - potentialLocationsTop(i).x);
            xDistances(xDistances>maxXDistance) = maxXDistance;
            unaries(j,:) = (maxXDistance - xDistances) / maxXDistance;
            unaries(j,:) = unaries(j,:) .* (potentialLocationsTop(i).y>minY)';
            
            % get pairwise potentials
            if i>1
                % get ind of last detection frame for object j
                if ~isnan(labels(i-1, j)) 
                    prevFrame = i-1;
                else
                    prevFrame = find(~isnan(labels(:,j)), 1, 'last');
                    wasOccluded(j) = 1;
                end

                % get label at previous dection frame
                prevLabel = labels(prevFrame, j);
                
                if isempty(prevFrame)
                    pairwise(j,:) = 0;
                else
                    pairwise(j,:) = getPairwisePotentials(potentialLocationsTop(i).x, potentialLocationsTop(i).y,...
                        potentialLocationsTop(prevFrame).x(prevLabel), potentialLocationsTop(prevFrame).y(prevLabel),...
                        frameTimes(i)-frameTimes(prevFrame), maxVel);
                    valid(j, pairwise(j,:)==0) = 0;
                end
            end
        end    
    end
    
    % get track scores (from svm)
    trackScores = repmat(potentialLocationsTop(i).scores', objectNum, 1);
    
    % get lowness of potential locations
    lowness = repmat(potentialLocationsTop(i).y' / vid.Height, objectNum, 1);

    % find best labels
    scores = unariesWeight.*unaries + pairwiseWeight.*pairwise + scoreWeight.*trackScores + lownessWeight.*lowness;
    scores = scores .* (unaries>0);
    scores = scores .* valid;
    labels(i,:) = getBestLabels(scores, objectNum, wasOccluded);

    % only keep labeled locations
    for j = 1:objectNum
        if isnan(labels(i,j))
            locationsTop.x(i,j) = nan;
            locationsTop.z(i,j) = nan;
        else
            locationsTop.x(i,j) = potentialLocationsTop(i).x(labels(i,j));
            locationsTop.z(i,j) = potentialLocationsTop(i).y(labels(i,j));
        end
    end
end
    
   
save([dataDir 'locationsTop.mat'], 'locationsTop');



% show tracking

if showTracking
    startFrame = 1;
    showPotentialLocations = true;
    showLines = false;
    
    if showLines
        showLocations(vid, potentialLocationsTop, labels, showPotentialLocations, .04, {[-1 -1], [-1 -1], [-1 -1], [-1 -1]}, startFrame, locationsBot);
    else
        showLocations(vid, potentialLocationsTop, labels, showPotentialLocations, .04, {[-1 -1], [-1 -1], [-1 -1], [-1 -1]}, startFrame);
    end
end





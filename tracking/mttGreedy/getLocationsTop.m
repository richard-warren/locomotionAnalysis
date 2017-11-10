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
maxXDistance = 15;
xOccludeBuffer = 40;
minY = 70;

unariesWeight = 1;
pairwiseWeight = 100;
scoreWeight = 0;


% initializations
labels = nan(length(potentialLocationsTop), objectNum);
locationsTop = struct();
vid = VideoReader(vidFile);
frameTimes = 0:.004:.004*(vid.NumberOfFrames-1); % temp, these are fake timestamps! oh shit!



% iterate through remaining frames

for i = 1:vid.NumberOfFrames
    
    unaries = zeros(objectNum, length(potentialLocationsTop(i).x));
    pairwise = zeros(objectNum, length(potentialLocationsTop(i).x));
    valid = ones(objectNum, length(potentialLocationsTop(i).x));
    wasOccluded = zeros(1, objectNum); % keeps track of whether the object was occluded in the previous frame (used in getBestLabels)
    
    
    for j = 1:objectNum
        
        % check if paw is occluded
        occludedByBins = abs(locationsBot(i).x(j) - locationsBot(i).x) < xOccludeBuffer & ...
                      (locationsBot(i).y > locationsBot(i).y(j));
        if any(occludedByBins)
            valid(j,:) = 0;
        
        % if not occluded, compute scores for all potential locations
        else
            
            % get unary potentials
            xDistances = abs(locationsBot(i).x(j) - potentialLocationsTop(i).x);
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

                % get label at last dection frame
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

    % find best labels
%     if i==23; keyboard; end
    scores = unariesWeight.*unaries + pairwiseWeight.*pairwise + scoreWeight.*trackScores;
    scores = scores .* (unaries>0);
    scores = scores .* valid;
    labels(i,:) = getBestLabels(scores, objectNum, wasOccluded);

    % only keep labeled locations
    for j = 1:objectNum
        if isnan(labels(i,j))
            locationsTop(i).x(j) = nan;
            locationsTop(i).y(j) = nan;
        else
            locationsTop(i).x(j) = potentialLocationsTop(i).x(labels(i,j));
            locationsTop(i).y(j) = potentialLocationsTop(i).y(labels(i,j));
        end
    end
end
    
   
save([dataDir 'locationsTop.mat'], 'locationsTop');





% show tracking
if showTracking
    startFrame = 1;
%     showLocations(vid, potentialLocationsTop, labels, .04, anchorPts, startFrame, locationsBot);
    showLocations(vid, potentialLocationsTop, labels, .04, {[-1 -1], [-1 -1], [-1 -1], [-1 -1]}, startFrame);
end

















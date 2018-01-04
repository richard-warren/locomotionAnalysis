function locationsTop = getLocationsTop(potentialLocationsTop, locationsBot, xLinearMapping, frameInds, obsPixPositions, frameTimeStamps, paws, fs)

% !!! need to document and make not shitty


% settings
objectNum = 4;
maxVel = 25 / .004;   % max velocity (pixels / sec)
maxXDistance = 25;    % max distance of x position in top from x position in bottom view
% xOccludeBuffer = 10;  % if paws in bottom view have x values within xOccludeBuffer of one another, the paw further away from the camera is treated as occluded
stanceHgt = 6;        % paws in stance (negative x velocity in bottom view) are assigned z values stanceHgt pixels above wheel (wheel defined by circRoiPts)
circRoiPts = [36 172; 224 122; 386 157];
stanceVel = -500;     % paws moving less than this vel are considered to be in stance (pixels / sec)
obsProximity = 60;    % if paw is within obsProximity pixels of obstacle, stance is no longer assumed when velocity is less than stanceVel

unariesWeight = 0;
pairwiseWeight = 0;
lownessWeight = 1;
scoreWeight = 0;


% initializations
labels = nan(length(potentialLocationsTop), objectNum);
locationsTop = struct();
[wheelRadius, wheelCenter] = fitCircle(circRoiPts - repmat([0 stanceHgt], 3, 1));

% fix x alignment for bottom view
locationsBot.x = locationsBot.x*xLinearMapping(1) + xLinearMapping(2);

% fix bottom view tracking (fixTracking fills short stretches of missing values and median filters)
locationsBot = fixTracking(locationsBot);

% get x velocities for bottom view tracking
locationsBot.xVel = nan(size(locationsBot.x));
for i = paws
    locationsBot.xVel(:,i) = getVelocity(locationsBot.x(:,i), .025, fs);
end



% iterate through all frameInds

for i = frameInds
    
    % report progress
    disp(i/length(locationsBot.x))
    
    % initializations
    unaries = zeros(objectNum, length(potentialLocationsTop(i).x));
    pairwise = zeros(objectNum, length(potentialLocationsTop(i).x));
    valid = ones(objectNum, length(potentialLocationsTop(i).x));
    wasOccluded = zeros(1, objectNum); % keeps track of whether the object was occluded in the previous frame (used in getBestLabels)
    
    
    for j = paws %1:objectNum
        
%         % check if paw is occluded
%         occludedByBins = abs(locationsBot.x(i,j) - locationsBot.x(i,:)) < xOccludeBuffer & ...
%                       (locationsBot.y(i,:) > locationsBot.y(i,j));
%         if any(occludedByBins)
%             valid(j,:) = 0;
%         
%         % if not occluded, compute scores for all potential locations
%         else
            
            % get unary potentials
            xDistances = abs(locationsBot.x(i,j) - potentialLocationsTop(i).x);
            xDistances(xDistances>maxXDistance) = maxXDistance;
            unaries(j,:) = (maxXDistance - xDistances) / maxXDistance;
            
            % get pairwise potentials
            if i>1
                % get ind of last detection frame for object j % !!! i think this is a bug... i cant just do i-1 because i-1 might not actually be a member of frameInds!!!
                if ~isnan(labels(i-1, j)) 
                    prevFrame = i-1;
                else
                    prevFrame = find(~isnan(labels(:,j)), 1, 'last');
                    wasOccluded(j) = 1;
                end

                % get label at previous dection frame
                prevLabel = labels(prevFrame, j);
                
                if isempty(prevFrame) || isempty(potentialLocationsTop(i).x)
                    pairwise(j,:) = 0;
                else
                    pairwise(j,:) = getPairwisePotentials(potentialLocationsTop(i).x, potentialLocationsTop(i).y,...
                        potentialLocationsTop(prevFrame).x(prevLabel), potentialLocationsTop(prevFrame).y(prevLabel),...
                        frameTimeStamps(i)-frameTimeStamps(prevFrame), maxVel);
                    valid(j, pairwise(j,:)==0) = 0;
                end
            end
%         end
    end
    
    % get track scores (from svm)
    trackScores = repmat(potentialLocationsTop(i).scores', objectNum, 1);
    
    % get lowness of potential locations
    lowness = repmat(potentialLocationsTop(i).y' / range(locationsBot.y(:)), objectNum, 1);

    % find best labels
    try
        scores = unariesWeight.*unaries + pairwiseWeight.*pairwise + scoreWeight.*trackScores + lownessWeight.*lowness;
        scores = scores .* (unaries>0);
        scores = scores .* valid;
    catch
        scores = []; % !!! this try/catch is a hack... need better way of handling frames where there are no potential locations
    end
    
    % !!! this if/then is temporary // it would be better if getBestLabels handled empty values itself // NEED TO LOOK INTO WHAT HAPPENS WHEN THERE ARE FEWER POTENTIAL OBJECTS THAN LOCATIONS... ARE NANS RETURNED?
    if isempty(scores)
        labels(i,:) = nan;
    else
        labels(i,:) = getBestLabels(scores, objectNum, wasOccluded);
    end
    

    % only keep labeled locations
    for j = paws %1:objectNum
        if isnan(labels(i,j))
            locationsTop.x(i,j) = nan;
            locationsTop.z(i,j) = nan;
        else
            locationsTop.x(i,j) = potentialLocationsTop(i).x(labels(i,j));
            locationsTop.z(i,j) = potentialLocationsTop(i).y(labels(i,j));
        end
        
        % if paw is in stance (negative velocity) make z location directly above the circle floor
        if locationsBot.xVel(i,j) < stanceVel && (abs(obsPixPositions(i) - locationsBot.x(i,j)) > obsProximity || isnan(obsPixPositions(i)))
            locationsTop.x(i,j) = locationsBot.x(i,j);
            locationsTop.z(i,j) = wheelCenter(2) - round(sqrt(wheelRadius^2 - (locationsBot.x(i,j)-wheelCenter(1))^2)); % !!! should replace this with pre-computed lookup table
        end
    end
end




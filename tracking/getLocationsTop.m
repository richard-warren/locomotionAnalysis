function locations = getLocationsTop(potentialLocationsTop, locationsBot, stanceBins, frameTimeStamps, wheelPoints)

% !!! need to document and make not shitty


% settings
objectNum = 4;
maxVel = 25 / .004;   % max velocity (pixels / sec)
maxXDistance = 30;    % max distance of x position in top from x position in bottom view
xOccludeBuffer = 25;  % if paws in bottom view have x values within xOccludeBuffer of one another, the paw further away from the camera is treated as occluded
stanceHgt = 5;        % paws in stance (negative x velocity in bottom view) are assigned z values stanceHgt pixels above wheel (wheel defined by circRoiPts)
maxHgt = 50;          % highness score are expressed as a fraction of maxHgt
minHgt = 6;          % only consider potential locations that are this many pixels above the circle (because we are only really looking for paws in swing bro)


unaryWeight = 1;
pairwiseWeight = 1;
highnessWeight = 1;


% initializations
locations.locationsRaw = nan(length(potentialLocationsTop), 2, objectNum);
[wheelRadius, wheelCenter] = fitCircle(wheelPoints);
wheelCenterOffset = wheelCenter - [0; stanceHgt];


% for all stanceBins frames locate the paw directly above the wheel at the x position recorded in locationsBot
for i = 1:4
    for j = find(stanceBins(:,i))'
        locations.locationsRaw(j,1,i) = locationsBot(j,1,i);
        locations.locationsRaw(j,2,i) = wheelCenterOffset(2) - round(sqrt(wheelRadius^2 - (locationsBot(j,1,i)-wheelCenterOffset(1))^2));
    end
end




% iterate through all frameInds
frameInds = find([potentialLocationsTop.isAnalyzed]);

for i = 1:length(frameInds)
    
    % sort paws from closest to furthest away from camera
    ys = squeeze(locationsBot(frameInds(i),2,:))';
    [~, pawSequence] = sort(ys, 'descend');
    if any(isnan(ys))
        nans = sum(isnan(ys));
        pawSequence = [pawSequence(1+nans:end) 1:nans];
    end
    alreadyTaken = false(length(potentialLocationsTop(frameInds(i)).x), 1);
    scores = nan(4, length(potentialLocationsTop(frameInds(i)).x));
    
    
    for j = 1:length(pawSequence)
        
        % check if paw is occluded
        isOccluded = any(abs(locationsBot(frameInds(i),1,pawSequence(j)) - locationsBot(frameInds(i),1,:)) < xOccludeBuffer &...
            (locationsBot(frameInds(i),2,:) > locationsBot(frameInds(i),2,pawSequence(j))));
        
        % only conduct analysis if paw is not occluded and if location has not been determined by stance analysis above
        if ~isOccluded && isnan(locations.locationsRaw(frameInds(i),2,pawSequence(j)))
            
            % get unary potentials
            xDistances = abs(locationsBot(frameInds(i),1,pawSequence(j)) - potentialLocationsTop(frameInds(i)).x);
            isTooFar = xDistances > maxXDistance;

            unaries = xDistances;
            unaries(isTooFar) = maxXDistance;
            unaries = (maxXDistance - unaries) / maxXDistance;
            unaries(isnan(unaries)) = 0;


            
            % get pairwise potentials

            % find last ind with detected paw
            lastInd = 1;
            for k = fliplr(1:i-1) % find last ind with detected frame
                if ~isnan(locations.locationsRaw(frameInds(k),1,pawSequence(j)))
                    lastInd = k;
                    break;
                end
            end

            dx = locations.locationsRaw(frameInds(lastInd),1,pawSequence(j)) - potentialLocationsTop(frameInds(i)).x;
            dz = locations.locationsRaw(frameInds(lastInd),2,pawSequence(j)) - potentialLocationsTop(frameInds(i)).z;
            dt = frameTimeStamps(frameInds(i)) - frameTimeStamps(frameInds(lastInd));
            vels = sqrt(dx.^2 + dz.^2) / dt;
            isTooFast = vels > maxVel;

            pairwise = vels;
            pairwise(isTooFast) = maxVel;
            pairwise = (maxVel - pairwise) / maxVel;
            pairwise(isnan(pairwise)) = 0;
            
            
            % get highness scores
            wheelZs = wheelCenter(2) - round(sqrt(wheelRadius^2 - (potentialLocationsTop(frameInds(i)).x-wheelCenter(1)).^2)); % !!! should replace this with pre-computed lookup table
            highness = wheelZs - potentialLocationsTop(frameInds(i)).z;
            tooLow = highness<minHgt;
            highness(highness>maxHgt) = maxHgt;
            highness = highness / maxHgt;

            
            % compute scores
            pawScores = (unaries * unaryWeight) + (pairwise * pairwiseWeight) + (highness * highnessWeight);
            pawScores(isTooFar | isTooFast | alreadyTaken | tooLow) = 0;
            scores(:, potentialLocationsTop(frameInds(i)).class==2) = 0; % remove potential locations not classified by cnn as a paw
            scores(pawSequence(j), :) = pawScores;
            

            [maxScore, maxInd] = max(scores(pawSequence(j), :));

            if maxScore>0
                locations.locationsRaw(frameInds(i),1,pawSequence(j)) = potentialLocationsTop(frameInds(i)).x(maxInd);
                locations.locationsRaw(frameInds(i),1,pawSequence(j)) = locationsBot(frameInds(i),1,pawSequence(j));
                locations.locationsRaw(frameInds(i),2,pawSequence(j)) = potentialLocationsTop(frameInds(i)).z(maxInd);
                alreadyTaken(maxInd) = true;
            end
        end
    end
end


% copy isAnalyzed and trialIdentities to locationsBot
locations.isAnalyzed = [potentialLocationsTop.isAnalyzed]';
locations.trialIdentities = [potentialLocationsTop.trialIdentities]';












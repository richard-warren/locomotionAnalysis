function locationsTop = getLocationsTop(potentialLocationsTop, locationsBot,...
    frameInds, wheelPositions, wheelTimes, wheelFs, mToPixFactor, obsPixPositions, frameTimeStamps, paws, fs)

% !!! need to document and make not shitty


% settings
objectNum = 4;
maxVel = 25 / .004;   % max velocity (pixels / sec)
maxXDistance = 30;    % max distance of x position in top from x position in bottom view
xOccludeBuffer = 25;  % if paws in bottom view have x values within xOccludeBuffer of one another, the paw further away from the camera is treated as occluded
stanceHgt = 6;        % paws in stance (negative x velocity in bottom view) are assigned z values stanceHgt pixels above wheel (wheel defined by circRoiPts)
circRoiPts = [36 155; 212 103; 374 132];
stanceVelDif = 1000;   % if paws paw is within this many pix/sec of wheel velocity (actually obs vel for now) then it is considered to be in stance IF length of this period exceeds stanceMin
stanceMin = .02;      % (s)
obsProximity = 60;    % if paw is within obsProximity pixels of obstacle, stance is no longer assumed when velocity is less than stanceVel
velTime = .04;       % amount of time to compute velocity over
maxHgt = 50;          % highness score are expressed as a fraction of maxHgt
minHgt = 6;          % only consider potential locations that are this many pixels above the circle (because we are only really looking for paws in swing bro)


unaryWeight = .5;
pairwiseWeight = .5;
highnessWeight = .5;


% initializations
locationsTop.x = nan(length(potentialLocationsTop), objectNum);
locationsTop.z = nan(length(potentialLocationsTop), objectNum);
[wheelRadius, wheelCenter] = fitCircle(circRoiPts);
wheelCenterOffset = wheelCenter - [0; stanceHgt]; % this guy is shifted up a little bit


% get x velocities for bottom view tracking
locationsBot.xVel = nan(size(locationsBot.x));
for i = paws
    locationsBot.xVel(:,i) = getVelocity(locationsBot.x(:,i), velTime, fs);
end


% get wheel velocity IN PIXELS
wheelVel = getVelocity(wheelPositions * mToPixFactor, velTime, wheelFs); % !!! need to change this to wheel speed...
wheelVel = interp1(wheelTimes, wheelVel, frameTimeStamps)';




% get paw locations for stance periods

for i = paws
    
    % !!! should really be doing this on a trial to trial basis // i dont think this code will handle well the transitions between trials...
    
    % get epoches where wheel vel and paw x vel are similar to one another
    matchedVelBins = abs(wheelVel(1:length(locationsBot.xVel(:,i))) - locationsBot.xVel(:,i)') < stanceVelDif;
    
    % exclude from stance consideration frames in which paw is close to obstacle
    nearObsBins = abs(obsPixPositions(1:length(locationsBot.x(:,i)))' - locationsBot.x(:,i)) < obsProximity;
    matchedVelBins(nearObsBins) = 0;
    
    startInds = find(diff(matchedVelBins) == 1) + 1;
    endInds = find(diff(matchedVelBins) == -1) + 1;

    % ensure that the first event is the beginning of an epoch and the last is the end of an epoch
    if endInds(1) < startInds(1); startInds = [1 startInds]; end
    if startInds(end) > endInds(end); endInds = [endInds length(matchedVelBins)]; end
    
    % only keep epochs that are long enough
    validStances = (frameTimeStamps(endInds) - frameTimeStamps(startInds)) > stanceMin;
    startInds = startInds(validStances);
    endInds = endInds(validStances);
    
    
    % store the coordinates of each paw during stance
    for j=1:length(startInds)
        for k = startInds(j):endInds(j)
            locationsTop.x(k,i) = locationsBot.x(k,i);
            locationsTop.z(k,i) = wheelCenterOffset(2) - round(sqrt(wheelRadius^2 - (locationsBot.x(k,i)-wheelCenterOffset(1))^2)); % !!! should replace this with pre-computed lookup table
        end
    end
end


% use the following to plot relationship between wheel and paw velocity // this can be used to determine stanceVelDif
% close all; figure;
% plot(wheelVel(1:length(locationsBot.xVel(:,i))));
% hold on; plot(locationsBot.xVel(:,i)');
% hold on; scatter(startInds, ones(1,length(startInds))*-1500)

% iterate through all frameInds

for i = 1:length(frameInds)
    
    % report progress
%     disp(i/length(frameInds))
    
    
    % sort paws from closest to furthest away from camera
    ys = locationsBot.y(frameInds(i),:);
    [~, pawSequence] = sort(ys, 'descend');
    if any(isnan(ys))
        nans = sum(isnan(ys));
        pawSequence = [pawSequence(1+nans:end) 1:nans];
    end
    alreadyTaken = false(length(potentialLocationsTop(frameInds(i)).x), 1);
    scores = nan(4, length(potentialLocationsTop(frameInds(i)).x));
    
    
    for j = 1:length(pawSequence)
        
        % check if paw is occluded
        isOccluded = any(abs(locationsBot.x(frameInds(i), pawSequence(j)) - locationsBot.x(frameInds(i),:)) < xOccludeBuffer &...
            (locationsBot.y(frameInds(i),:) > locationsBot.y(frameInds(i), pawSequence(j))));
        
        
        if ~isOccluded && isnan(locationsTop.z(frameInds(i), pawSequence(j))) % only conduct analysis if location has not been determined by stance analysis above
            
            % get unary potentials
            xDistances = abs(locationsBot.x(frameInds(i), pawSequence(j)) - potentialLocationsTop(frameInds(i)).x);
            isTooFar = xDistances > maxXDistance;

            unaries = xDistances;
            unaries(isTooFar) = maxXDistance;
            unaries = (maxXDistance - unaries) / maxXDistance;
            unaries(isnan(unaries)) = 0;


            
            % get pairwise potentials

            % find last ind with detected paw
            lastInd = 1;
            for k = fliplr(1:i-1) % find last ind with detected frame
                if ~isnan(locationsTop.x(frameInds(k), pawSequence(j)))
                    lastInd = k;
                    break;
                end
            end

            dx = locationsTop.x(frameInds(lastInd), pawSequence(j)) - potentialLocationsTop(frameInds(i)).x;
            dz = locationsTop.z(frameInds(lastInd), pawSequence(j)) - potentialLocationsTop(frameInds(i)).z;
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
%             if frameInds(i) == 22188; keyboard; end


            
            % compute scores
            pawScores = (unaries * unaryWeight) + (pairwise * pairwiseWeight) + (highness * highnessWeight);
            pawScores(isTooFar | isTooFast | alreadyTaken | tooLow) = 0;
            scores(pawSequence(j), :) = pawScores;

            [maxScore, maxInd] = max(scores(pawSequence(j), :));

            if maxScore>0
                locationsTop.x(frameInds(i), pawSequence(j)) = potentialLocationsTop(frameInds(i)).x(maxInd);
                locationsTop.z(frameInds(i), pawSequence(j)) = potentialLocationsTop(frameInds(i)).z(maxInd);
                alreadyTaken(maxInd) = true;
            end
        end
    end
end













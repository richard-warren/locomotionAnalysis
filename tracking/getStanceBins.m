function stanceBins = getStanceBins(vidTop, wheelPoints, xLocations, trialIdentities, fs, mToPixFactor, ...
    wheelPositions, wheelTimes, targetFs, frameTimeStamps)

% !!! needs documentation, but generally returns n x 4 binary vector recording whether 4 paws are in stance
% does this by comparing wheel and paw velocities // when they are matched for > stanceMin seconds, we call it stance

% settings
showAnalysis = false;
stanceVelDif = 1000;   % if paws paw is within this many pix/sec of wheel velocity (actually obs vel for now) then it is considered to be in stance IF length of this period exceeds stanceMin
stanceMin = .02;       % (s)
velTime = .02;         % amount of time to compute velocity over

% settings for checking that foot is actually touching flor
subFrameUpDown = [10 6];
subFrameLeftRight = [12 12];
stanceThresh = 80;

% get wheel parameters
[wheelRadius, wheelCenter] = fitCircle(wheelPoints);

% get x velocities for bottom view tracking
xVel = nan(size(xLocations));
for i = 1:4
    xVel(:,i) = getVelocity(xLocations(:,i), velTime, fs);
end


% get wheel velocity IN PIXELS
wheelVel = getVelocity(wheelPositions * mToPixFactor, velTime, targetFs);
wheelVel = interp1(wheelTimes, wheelVel, frameTimeStamps)';


% get stance bins for each trial
stanceBinsUncorrected = false(size(xLocations));
stanceBins = false(size(xLocations));
isTouching = false(size(xLocations));
startEndInds = cell(1,4);
% !!! temp hack fix // wont need this line after potentialLocationsTop are reanalyzed st length of structure is same as number of frames in video
trialIdentities(end+1:length(frameTimeStamps)) = nan;



for i = unique(trialIdentities(~isnan(trialIdentities)))
    for j = 1:4

        % get epoches where wheel vel and paw x vel are similar to one another
        matchedVelBins = (abs(wheelVel - xVel(:,j)') < stanceVelDif) & trialIdentities==i;

        % exclude from stance consideration frames in which paw is close to obstacle
        % !!! need to check that excluding these lines of code doesn't create false stances when he is butting up against the obstacle
%         nearObsBins = abs(obsPixPositions' - locationsBot(:,1,j)) < obsProximity;
%         matchedVelBins(nearObsBins) = 0;

        startInds = find(diff(matchedVelBins) == 1) + 1;
        endInds = find(diff(matchedVelBins) == -1) + 1;

        % ensure that the first event is the beginning of an epoch and the last is the end of an epoch
        if endInds(1) < startInds(1); startInds = [1 startInds]; end
        if startInds(end) > endInds(end); endInds = [endInds length(matchedVelBins)]; end

        % only keep epochs that are long enough
        validStances = (frameTimeStamps(endInds) - frameTimeStamps(startInds)) > stanceMin;
        startInds = startInds(validStances);
        endInds = endInds(validStances);

        % store results
        for k = 1:length(startInds)
            stanceBinsUncorrected(startInds(k):endInds(k),j) = true;
        end
        startEndInds{j} = cat(2, startEndInds{j}, [startInds; endInds]);
        
    end
end



% prepare figure is showAnalysis
if showAnalysis
    figure; colormap gray; pimpFig
    
    % raw frame preview, with rectangle showing subframe region
    subplot(2,2,1:2)
    framePreview = imshow(randi(255, vidTop.Height, vidTop.Width));
    set(gca, 'visible', 'off', 'CLim', [0 255]);
    rect = rectangle(gca, 'Position', [0 0 0 0], 'EdgeColor', 'red');
    
    % raw subframe preview
    subplot(2,2,3)
    rawPreview = imshow(randi(255, sum(subFrameUpDown), sum(subFrameLeftRight)));
    set(gca, 'visible', 'off', 'CLim', [0 255]);
    
    % threshed subframe preview
    subplot(2,2,4)
    threshedPreview = imshow(false(diff(subFrameUpDown), diff(subFrameLeftRight)));
end




% for each detected stance, determine whether foot is touching wheel
for i = find(sum(stanceBinsUncorrected,2))'
    disp(i)
    
    frame = rgb2gray(read(vidTop, i));
    stancePaws = find(stanceBinsUncorrected(i,:));
    
    for j = stancePaws
        
        % ensure that foot is touching floor
        pawX = round(xLocations(i,j));
        wheelY = round(wheelCenter(2) - round(sqrt(wheelRadius^2 - (xLocations(i,j)-wheelCenter(1))^2)));
        subFrame = frame(max(wheelY-subFrameUpDown(1), 1) : min(wheelY+subFrameUpDown(2), size(frame,1)), ...
                         max(pawX-subFrameLeftRight(1), 1) : min(pawX+subFrameLeftRight(2), size(frame,2)));
        
        regions = bwlabel(subFrame<stanceThresh);
        leftSideRegions = unique(regions(regions(:,1)~=0,1));
        rightSideRegions = unique(regions(regions(:,end)~=0,end));
        isTouching(i,j) = ~any(intersect(leftSideRegions, rightSideRegions));
        
        
        if showAnalysis
            set(framePreview, 'CData', frame);
            set(rawPreview, 'CData', subFrame)
            set(threshedPreview, 'CData', ~(subFrame<stanceThresh))
            set(rect, 'Position', [pawX-subFrameLeftRight(1), wheelY-subFrameUpDown(2), ...
                sum(subFrameLeftRight), sum(subFrameUpDown)]);
            disp(i)
            waitforbuttonpress
        end
    end
end


% fine tune stance start and end inds
% (go through all atance start and end inds, and adjust them s.t. each stance starts with the first isTouching frame in the stance and ends with last isTouching frame)
for i = 1:4
    for j = 1:length(startEndInds{i})
        
        currentStanceBins = 1:size(isTouching,1)>=startEndInds{i}(1,j) & 1:size(isTouching,1)<=startEndInds{i}(2,j);
        firstTouchingInd = find(isTouching(:,i)' & currentStanceBins, 1, 'first'); % first ind within uncorrected stance in which paw is touching
        lastTouchingInd = find(isTouching(:,i)' & currentStanceBins, 1, 'last'); % last ind within uncorrected stance in which paw is touching
        
        stanceBins(firstTouchingInd:lastTouchingInd, i) = true;
    end
end






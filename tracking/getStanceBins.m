function stanceBins = getStanceBins(vidTop, wheelPoints, xLocations, trialIdentities, fs, mToPixFactor, ...
    wheelPositions, wheelTimes, targetFs, frameTimeStamps)

% !!! needs documentation, but generally returns n x 4 binary vector recording whether 4 paws are in stance
% does this by comparing wheel and paw velocities // when they are matched for > stanceMin seconds, we call it stance

% settings
showAnalysis = false;
stanceVelDif = 1000;   % if paws paw is within this many pix/sec of wheel velocity (actually obs vel for now) then it is considered to be in stance IF length of this period exceeds stanceMin
stanceSwingMin = .02;       % (s) min time for either stance or swing
velTime = .02;         % amount of time to compute velocity over

% settings for checking that foot is actually touching flor
subFrameUpDown = [12 6];
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



% !!! add way to ignore excluded trials


for i = unique(trialIdentities(~isnan(trialIdentities)))
    for j = 1:4

        % get epoches where wheel vel and paw x vel are similar to one another
        matchedVelBins = (abs(wheelVel - xVel(:,j)') < stanceVelDif) & trialIdentities==i;

        % exclude from stance consideration frames in which paw is close to obstacle
        % !!! need to check that excluding these lines of code doesn't create false stances when he is butting up against the obstacle
        % !!! could change this so the subFrame computed below shift to the right or left if the paw is ahead of or behaind the obs, s.t. obs never sppears in subFrame...
%         nearObsBins = abs(obsPixPositions' - locationsBot(:,1,j)) < obsProximity;
%         matchedVelBins(nearObsBins) = 0;
        

        stanceBinsUncorrected(matchedVelBins,j) = true;
        switchInds = find(diff(matchedVelBins)~=0)+1;
        
        % remove short swing and stances
        for k = 1:(length(switchInds)-1)
            if ~isnan(switchInds(k)) % if switch ind hasn't been removed in previous iteration of the loop
                if (frameTimeStamps(switchInds(k+1)) - frameTimeStamps(switchInds(k))) < stanceSwingMin
                    inds = switchInds(k):(switchInds(k+1)-1);
                    stanceBinsUncorrected(inds,j) = ~stanceBinsUncorrected(inds,j);
                    switchInds(k+1) = nan;
                end
            end
        end
    end
end



% prepare figure is showAnalysis
if showAnalysis
    figure; colormap gray; pimpFig; set(gcf, 'color', [.5 .5 .5])
    
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
w = waitbar(0, 'getting stance bins...', 'position', [1500 50 270 56.2500]);

for i = find(sum(stanceBinsUncorrected,2))'
    
    frame = rgb2gray(read(vidTop, i));
    stancePaws = find(stanceBinsUncorrected(i,:));
    
    for j = stancePaws
        
        % ensure that foot is touching floor
        pawX = round(xLocations(i,j));
        wheelY = round(wheelCenter(2) - round(sqrt(wheelRadius^2 - (xLocations(i,j)-wheelCenter(1))^2)));
        subFrame = frame(max(wheelY-subFrameUpDown(1), 1) : min(wheelY+subFrameUpDown(2), size(frame,1)), ...
                         max(pawX-subFrameLeftRight(1), 1) : min(pawX+subFrameLeftRight(2), size(frame,2)));
        
        if ~isempty(subFrame)
            regions = bwlabel(subFrame<stanceThresh);
            leftSideRegions = unique(regions(regions(:,1)~=0,1));
            rightSideRegions = unique(regions(regions(:,end)~=0,end));
            isTouching(i,j) = ~any(intersect(leftSideRegions, rightSideRegions));
        else
            fprintf('  problem with trial %i\n', i)
            keyboard
        end
        
        
        if showAnalysis && (i==32955) %  || i==23574
            set(framePreview, 'CData', frame);
            set(rawPreview, 'CData', subFrame)
            set(threshedPreview, 'CData', ~(subFrame<stanceThresh))
            set(rect, 'Position', [pawX-subFrameLeftRight(1), wheelY-subFrameUpDown(2), ...
                sum(subFrameLeftRight), sum(subFrameUpDown)]);
            disp(i)
            waitforbuttonpress
        end
    end
    
    waitbar(i/size(stanceBins,1))
end



% fine tune stance start and end inds
% (go through all atance start and end inds, and adjust them s.t. each stance starts with the first isTouching frame in the stance and ends with last isTouching frame)
for i = 1:4
    
    startInds = find(diff(stanceBinsUncorrected(:,i))==1) + 1;
    endInds = find(diff(stanceBinsUncorrected(:,i))==-1);
    if endInds(1)<startInds(1); endInds = endInds(2:end); end
    if startInds(end)>endInds(end); startInds = startInds(1:end-2); end
    
    for j = 1:length(startInds)
        
        currentStanceBins = 1:size(isTouching,1)>=startInds(j) & 1:size(isTouching,1)<=endInds(j);
        firstTouchingInd = find(isTouching(:,i)' & currentStanceBins, 1, 'first'); % first ind within uncorrected stance in which paw is touching
        lastTouchingInd = find(isTouching(:,i)' & currentStanceBins, 1, 'last'); % last ind within uncorrected stance in which paw is touching
        
        stanceBins(firstTouchingInd:lastTouchingInd, i) = true;
    end
end

close(w)



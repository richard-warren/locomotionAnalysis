function stanceBins = getStanceBins(frameTimeStamps, locationsTopPaws, wheelPositions, wheelTimes, wheelCenter, wheelRadius, vidFs, mToPixFactor)



% settings
stanceVelDif = 1000;   % if paws paw is within this many pix/sec of wheel velocity (actually obs vel for now) then it is considered to be in stance IF length of this period exceeds stanceMin
stanceSwingMin = .02;       % (s) min time for either stance or swing
velTime = .02;         % amount of time to compute velocity over
stanceMaxHgt = 12;



% get x velocities for bottom view tracking
xVel = nan(size(locationsTopPaws,1),4);
for i = 1:4
    xVel(:,i) = getVelocity(locationsTopPaws(:,1,i), velTime, vidFs);
end

% get wheel velocity IN PIXELS
wheelVel = getVelocity(wheelPositions * mToPixFactor, velTime, 1/(median(diff(wheelTimes))));
wheelVel = interp1(wheelTimes, wheelVel, frameTimeStamps)';






% determine when paw is close to wheel
xLocations = squeeze(locationsTopPaws(:,1,:));
yMaxes = wheelCenter(2) - round(sqrt(wheelRadius^2 - (xLocations-wheelCenter(1)).^2));
isCloseToWheel = squeeze(locationsTopPaws(:,2,:)) > (yMaxes-stanceMaxHgt); % remember that y values are flipped upside down
% pawHeights = squeeze(yMaxes - locationsTopPaws(:,2,:));
% figure; histogram(pawHeights(stanceBinsUncorrected)); hold on; histogram(pawHeights(~stanceBinsUncorrected));




% determine when paw is in stance
stanceBins = false(size(locationsTopPaws,1),4);


for i = 1:4
    % get epoches where wheel vel and paw x vel are similar to one another
    potentialStanceBins = (abs(wheelVel - xVel(:,i)') < stanceVelDif) & isCloseToWheel(:,i)';

    stanceBins(potentialStanceBins,i) = true;
    switchInds = find(diff(potentialStanceBins)~=0)+1; % when paw switches from swing to stance or stance to swing

    % remove short swing and stances
    for j = 1:(length(switchInds)-1)
        if ~isnan(switchInds(j)) % if switch ind hasn't been removed in previous iteration of the loop
            if (frameTimeStamps(switchInds(j+1)) - frameTimeStamps(switchInds(j))) < stanceSwingMin
                inds = switchInds(j):(switchInds(j+1)-1);
                stanceBins(inds,i) = ~stanceBins(inds,i);
                switchInds(j+1) = nan;
            end
        end
    end
end


    




% % fine tune stance start and end inds
% % (go through all atance start and end inds, and adjust them s.t. each stance starts with the first isCloseToWheel frame in the stance and ends with last isTouching frame)
% for i = 1:4
%     
%     startInds = find(diff(stanceBinsUncorrected(:,i))==1) + 1;
%     endInds = find(diff(stanceBinsUncorrected(:,i))==-1);
%     if endInds(1)<startInds(1); endInds = endInds(2:end); end
%     if startInds(end)>endInds(end); startInds = startInds(1:end-2); end
%     
%     for j = 1:length(startInds)
%         
%         currentStanceBins = 1:size(isCloseToWheel,1)>=startInds(j) & 1:size(isCloseToWheel,1)<=endInds(j);
%         firstTouchingInd = find(isCloseToWheel(:,i)' & currentStanceBins, 1, 'first'); % first ind within uncorrected stance in which paw is touching
%         lastTouchingInd = find(isCloseToWheel(:,i)' & currentStanceBins, 1, 'last'); % last ind within uncorrected stance in which paw is touching
%         
%         stanceBins(firstTouchingInd:lastTouchingInd, i) = true;
%     end
% end





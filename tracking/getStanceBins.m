function stanceBins = getStanceBins(xLocations, trialIdentities, fs, mToPixFactor, ...
    wheelPositions, wheelTimes, targetFs, frameTimeStamps)

% !!! needs documentation, but generally returns n x 4 binary vector recording whether 4 paws are in stance
% does this by comparing wheel and paw velocities // when they are matched for > stanceMin seconds, we call it stance

% settings
stanceVelDif = 1000;   % if paws paw is within this many pix/sec of wheel velocity (actually obs vel for now) then it is considered to be in stance IF length of this period exceeds stanceMin
stanceMin = .02;       % (s)
velTime = .02;         % amount of time to compute velocity over


% get x velocities for bottom view tracking
xVel = nan(size(xLocations));
for i = 1:4
    xVel(:,i) = getVelocity(xLocations(:,i), velTime, fs);
end


% get wheel velocity IN PIXELS
wheelVel = getVelocity(wheelPositions * mToPixFactor, velTime, targetFs);
wheelVel = interp1(wheelTimes, wheelVel, frameTimeStamps)';


% get stance bins for each trial
stanceBins = false(size(xLocations));

% !!! temp hack fix // wont need this line after potentialLocationsTop are reanalyzed st length of structure is same as number of frames in video
trialIdentities(end+1:length(frameTimeStamps)) = nan;


for i = unique(trialIdentities(~isnan(trialIdentities)))
    
    for j = 1:4

        % get epoches where wheel vel and paw x vel are similar to one another
        try
        matchedVelBins = (abs(wheelVel - xVel(:,j)') < stanceVelDif) & ...
                          trialIdentities==i;
        catch; keyboard; end

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
            stanceBins(startInds(k):endInds(k),j) = true;
        end
        
    end
end


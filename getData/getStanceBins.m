function stanceBins = getStanceBins(...
    frameTimeStamps, locationsTopPaws, wheelPositions, wheelTimes, wheelCenter, wheelRadius, vidFs, pixelsPerM)

% returns a (frameNum X 4) binary matrix recording whether each paw is in
% stance for each frame // does this by ensuring the paw is within
% stanceMaxHgt of the wheel AND the x velocity is within stanceVelDif m/s
% of the wheel velocity // it then applies median filtering to debounce

% settings
stanceVelDif = .2;       % (m/s) if paws paw is within this many pix/sec of wheel velocity then it is considered to be in stance IF length of this period exceeds stanceMin
stanceMaxHgt = .005;     % (m) max vertical distance from wheel for paw to be considered potentially touching wheel
medFiltDuration = .02;   % (s) length of temporal filtering kernel
velTime = .03;           % amount of time to compute velocity over


% get paw and wheel velocities
xVel = nan(size(locationsTopPaws,1),4);
for i = 1:4; xVel(:,i) = getVelocity(locationsTopPaws(:,1,i), velTime, vidFs); end
xVel = xVel / pixelsPerM;

wheelPosInterp = interp1(wheelTimes, -wheelPositions, frameTimeStamps);  % negative wheelPositions because wheel is rotating towards the left of the screen
wheelVel = getVelocity(wheelPosInterp, velTime, vidFs);

% determine paw and wheel velocities are close to eachother
isVelClose = abs(xVel - wheelVel')<stanceVelDif;

% determine when paw is close to wheel
xLocations = squeeze(locationsTopPaws(:,1,:));
yMaxes = wheelCenter(2) - round(sqrt(wheelRadius^2 - (xLocations-wheelCenter(1)).^2));
isPawClose = squeeze(locationsTopPaws(:,2,:)) > (yMaxes-(stanceMaxHgt*pixelsPerM));  % remember that y values are flipped upside down


% determine when paw is in stance
stanceBins = isVelClose & isPawClose;
stanceBins = medfilt2(stanceBins, [round(medFiltDuration*vidFs) 1]);






% get x velocities for bottom view tracking
locationsBot.xVel = nan(size(locationsBot.x));
locationsBotFixed = fixTracking(locationsBot);

for i = paws
    locationsBot.xVel(:,i) = getVelocity(locationsBotFixed.x(:,i), .025, fs);
end

% get wheel velocity
wheelVel = getVelocity(obsPixPositions, .025, fs);


%% plot paw and wheel velocities
close all; figure; plot(locationsBot.xVel(:,1)); hold on; plot(wheelVel); pimpFig
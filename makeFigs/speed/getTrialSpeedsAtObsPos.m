function getTrialSpeedsAtObsPos(obsPos, wheelPositions, obsPositions, obsTimes, obsOnTimes, speedTime)

% computes wheel velocity when obstacle is at obsPos


% initializations
vel = getVelocity(wheelPositions, speedTime, targetFs);
trialVels = nan(1, length(obsOnTimes));

for i = 1:length(obsOnTimes)
    
    % get trial vel and times
    obsAtPosTime = obsTimes(find(obsTimes>=obsOnTimes(i) & obsPositions>=obsPos, 1, 'first'));
    trialVels(i) = interp1(wheelTimes, vel, obsAtPosTime);
    
end
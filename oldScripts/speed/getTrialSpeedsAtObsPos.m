function sessionVels = getTrialSpeedsAtObsPos(obsPos, wheelPositions, wheelTimes, obsPositions, obsTimes, obsOnTimes, speedTime, targetFs)

% computes wheel velocity when obstacle is at obsPos


% initializations
vel = getVelocity(wheelPositions, speedTime, targetFs);
sessionVels = nan(1, length(obsOnTimes));

for i = 1:length(obsOnTimes)
    
    % get trial vel and times
    obsAtPosTime = obsTimes(find(obsTimes>=obsOnTimes(i) & obsPositions>=obsPos, 1, 'first'));
    sessionVels(i) = interp1(wheelTimes, vel, obsAtPosTime);
    
end
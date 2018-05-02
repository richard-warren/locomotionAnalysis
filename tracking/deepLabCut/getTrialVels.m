function trialVels = getTrialVels(velPrePost, obsOnTimes, obsTimes, obsPositions)

% this code assumes that obsPositions have already been normalized s.t. 0 is where obs is beneathe animal's nose

trialVels = nan(1, length(obsOnTimes));

for i = 1:length(obsOnTimes)
    
    % get trial velocity
    startInd = find(obsTimes>obsOnTimes(i) & obsPositions>-velPrePost(1), 1, 'first');
    endInd = find(obsTimes>obsOnTimes(i) & obsPositions>velPrePost(2), 1, 'first');
    trialVels(i) = (obsPositions(endInd) - obsPositions(startInd)) / (obsTimes(endInd) - obsTimes(startInd));
    
end
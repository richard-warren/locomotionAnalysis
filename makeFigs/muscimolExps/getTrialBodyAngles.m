function trialBodyAngles = getTrialBodyAngles(bodyAngles, obsOnTimes, obsOffTimes, frameTimeStamps, breaks)

% computes the average body angle for each trial while obs is on, ommitting
% frames that occur after a wheel break

trialBodyAngles = nan(1,length(obsOnTimes));

for i = 1:length(obsOnTimes)
    
    trialInds = frameTimeStamps>obsOnTimes(i) & frameTimeStamps<obsOffTimes(i);
    if any(breaks.times>obsOnTimes(i) & breaks.times<obsOffTimes(i))
        breakTime = breaks.times(breaks.times>obsOnTimes(i) & breaks.times<obsOffTimes(i));
        try
        trialInds = trialInds(frameTimeStamps(trialInds)<breakTime(1));
        catch; keyboard; end
    end
    
    trialBodyAngles(i) = nanmedian(bodyAngles(trialInds));
end
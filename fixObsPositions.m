function obsPositionsFixed = fixObsPositions(obsPositions, obsTimes, obsOnTimes)

% fixes the drift in the rotary encoder reading position of obstacle stepper motor
% this drift could in principle result from either from lost ticks in the endoer OR belt skidding, but it is likely the former
% for every obstacle trial, it finds the minimum position (which is reached when the platform re-zeros to the end stop after each trial)
% it then subtracts this position for all other positions, effectively settings to 0 the lowest position in each trial
%
% input      obsPositions:        raw obstacle position readings (m), read from rotary decoder coupled to stepper motor shaft (actually the idler on the opposite side of the track)
%            obsTimes:            timestamps for obsPositions
%            obsOnTimes:          the times at which the obstacle turns on
%
% output     obsPositionsFixed:   obsPositions with the trial to trial signal drift removed

obsPositionsFixed = obsPositions;

for i = 1:length(obsOnTimes)
    
    if i<length(obsOnTimes)
        endTime = obsOnTimes(i+1);
    else
        endTime = max(obsTimes);
    end
    
    trialInds = (obsTimes>obsOnTimes(i)) & (obsTimes<=endTime);
    minPos = min(obsPositions(trialInds));
    obsPositionsFixed(trialInds) = obsPositionsFixed(trialInds) - minPos;
    
end
function [frameInds, trialIdentities, trialVels] = getTrialFrameInds(minVel, obsPrePost, velPrePost, frameTimeStamps,...
    wheelPositions, wheelTimes, obsPositions, obsTimes, obsOnTimes)

% returns vector containing frameInds for all trials // only includes trials where mouse is running faster than minVel
% minVel is calculated between velPositions, which should be the start and end positions of the mouse running over obs
% frameInds include all frameInds while obstacle is on, but also obsPrePost(1) meters before the obs turns on
% and obsprePost(2) meters after the obs turns off // setting this to [0 0] means only inds with obs on are included
% also saves trialIdentities, which for each frameInd records the trial number

% !!! now obsPrePost actually determines how many m before and after obs is at tip of nose
% this code assumes that obsPositions have already been normalized s.t. 0 is where obs is beneathe animal's nose

% initializations
frameInds = [];
trialIdentities = [];
trialVels = nan(1, length(obsOnTimes));


for i = 1:length(obsOnTimes)
    
    % get trial velocity
    startInd = find(obsTimes>obsOnTimes(i) & obsPositions>-velPrePost(1), 1, 'first');
    endInd = find(obsTimes>obsOnTimes(i) & obsPositions>velPrePost(2), 1, 'first');
    trialVels(i) = (obsPositions(endInd) - obsPositions(startInd)) / (obsTimes(endInd) - obsTimes(startInd));
    
    if trialVels(i) > minVel
    
        % get positions of wheel at moment obs is at wheel center
        noseTime = obsTimes(find(obsTimes>obsOnTimes(i) & obsPositions>=0, 1, 'first'));
        wheelPos = interp1(wheelTimes, wheelPositions, noseTime);

        % add and subtract obsPrePost to these values
        startWheelPos = wheelPos - obsPrePost(1);
        endWheelPos = wheelPos + obsPrePost(2);

        % find times corresponding to these start and stop positions
        startTime = wheelTimes(find(wheelPositions>=startWheelPos, 1, 'first'));
        endTime = wheelTimes(find(wheelPositions>=endWheelPos, 1, 'first'));

        % getFrameInds
        newFrames = find(frameTimeStamps>=startTime & frameTimeStamps<=endTime)';
        frameInds = [frameInds newFrames];
        
        % store trial identities
        trialIdentities = [trialIdentities ones(1,length(newFrames))*i];
        
    end
end








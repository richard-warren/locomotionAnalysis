function [frameInds, trialIdentities, trialVels] = getTrialFrameInds(minVel, obsPrePost, velPositions, frameTimeStamps,...
    wheelPositions, wheelTimes, obsPositions, obsTimes, obsOnTimes, obsOffTimes)

% returns vector containing frameInds for all trials // only includes trials where mouse is running faster than minVel
% minVel is calculated between velPositions, which should be the start and end positions of the mouse running over obs
% frameInds include all frameInds while obstacle is on, but also obsPrePost(1) meters before the obs turns on
% and obsprePost(1) meters after the obs turns off // setting this to [0 0] means only inds with obs on are included
% also saves trialIdentities, which for each frameInd records the trial number

% temp
% session = '180122_001';
% load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], 'wheelPositions',...
%     'wheelTimes', 'obsOnTimes', 'obsOffTimes', 'frameTimeStamps', 'obsPositions', 'obsTimes')
% obsPositions = fixObsPositions(obsPositions, obsTimes, obsOnTimes);
% minVel = .4;
% obsPrePost = [.2 .2];
% velPositions = [-.08 .08] + 0.3820;

frameInds = [];
trialIdentities = [];
trialVels = nan(1, length(obsOnTimes));


for i = 1:length(obsOnTimes)
    
    % get trial velocity
    startInd = find(obsTimes>obsOnTimes(i) & obsPositions>velPositions(1), 1, 'first');
    endInd = find(obsTimes>obsOnTimes(i) & obsPositions>velPositions(2), 1, 'first');
    trialVels(i) = (obsPositions(endInd) - obsPositions(startInd)) / (obsTimes(endInd) - obsTimes(startInd));
    
    if trialVels(i) > minVel
    
        % get positions of wheel at moments obs turn on and off
        startWheelPos = interp1(wheelTimes, wheelPositions, obsOnTimes(i));
        endWheelPos = interp1(wheelTimes, wheelPositions, obsOffTimes(i));

        % add and subtract obsPrePost to these values
        startWheelPos = startWheelPos - obsPrePost(1);
        endWheelPos = endWheelPos + obsPrePost(2);

        % find times corresponding to these start and stop positions
        startTime = wheelTimes(find(wheelPositions>=startWheelPos, 1, 'first'));
        endTime = wheelTimes(find(wheelPositions>=endWheelPos, 1, 'first'));

        % getFrameInds
        frameInds = [frameInds find(frameTimeStamps>=startTime & frameTimeStamps<=endTime)'];
        
        % store trial identities
        trialIdentities = [trialIdentities ones(1,length(frameInds))*i];
        
    end
    
end
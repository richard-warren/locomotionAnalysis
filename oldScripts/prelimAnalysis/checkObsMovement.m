% function checkObsMovement

% a principled way to do this: if the obs light is on, check pos tracking only while it is ACTUALLY on...
% if light is not on, only check pos tracking within some positional range around the wheel (because he cant touch it elsewhere and it more or less doesnt matter)
% also, only check up to moment of wheel break... its to be expected that the tracking will fail at these times

load('C:\Users\Rick\Google Drive\columbia\obstacleData\sessions\171019_002\runAnalyzed.mat')
obsPositions = fixObsPositions(obsPositions, obsTimes, obsOnTimes);

velTime = .1;
obsVel = getVelocity(obsPositions, velTime, targetFs);
wheelVel = getVelocity(wheelPositions, velTime, targetFs);
%%


close all; figure;
cmap = winter(length(obsOnTimes));

for i = 1:length(obsOnTimes)

    % get obs and wheel indices
    obsInds = find(obsTimes>obsOnTimes(i) & obsTimes<obsOffTimes(i));
    wheelInds = find(wheelTimes>obsOnTimes(i) & wheelTimes<obsOffTimes(i));
    
    % ensure obs and wheel ind vectors are same length
    trialLength = min(length(obsInds), length(wheelInds));
    obsInds = obsInds(1:trialLength);
    wheelInds = wheelInds(1:trialLength);
    
    
    % get trial postiions (for obs and wheel)
    obsTrial = obsPositions(obsInds) - obsPositions(obsInds(end));
    wheelTrial = wheelPositions(wheelInds) - wheelPositions(wheelInds(end));
    
    % get trial velocities (for obs and wheel)
    obsTrialVel = obsVel(obsInds);
    wheelTrialVel = wheelVel(wheelInds);

    % get position and velocity differences
    trialPosDiff = (obsTrial - wheelTrial);
    trialVelDiff = (obsTrialVel - wheelTrialVel);

    % make time vector
    times = linspace(0, obsOffTimes(i)-obsOnTimes(i), trialLength);
    
    % plot trial
    subplot(2,1,1)
    plot(times, trialPosDiff, 'color', cmap(i,:)); hold on
    
    subplot(2,1,2)
    plot(times, trialVelDiff, 'color', cmap(i,:)); hold on

end


% set(gca, 'xlim', [0 median(obsOffTimes - obsOnTimes)])
for i=1:2
    subplot(2,1,i)
    set(gca, 'xlim', [0 median(obsOffTimes - obsOnTimes)])
end
pimpFig
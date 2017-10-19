% function checkObsMovement

obsPositions = fixObsPositions(obsPositions, obsTimes, obsOnTimes);
%%

figure;
cmap = winter(length(obsOnTimes));
for i = 1:length(obsOnTimes)



    obsTrial = obsPositions(obsTimes>obsOnTimes(i) & obsTimes<obsOffTimes(i));
    wheelTrial = wheelPositions(wheelTimes>obsOnTimes(i) & wheelTimes<obsOffTimes(i));

    trialLength = min(length(obsTrial), length(wheelTrial));
    obsTrial = obsTrial(1:trialLength) - obsTrial(end);
    wheelTrial = wheelTrial(1:trialLength) - wheelTrial(end);

    times = linspace(0, obsOffTimes(i)-obsOnTimes(i), trialLength);
    trialDiff = abs(obsTrial - wheelTrial);

    plot(times, trialDiff, 'color', cmap(i,:)); hold on

end

pimpFig
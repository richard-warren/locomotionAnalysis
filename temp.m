%% measure times of trials and other little thangs

fprintf('loading... '); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'baseline_data.mat'), 'data'); disp('baseline data loaded!')

flat = flattenData(data, {'wiskContactTimes', 'obsOnTimes', 'obsOffTimes', 'lightOnTimes', 'isLightOn'});
flat = flat(~[flat.isLightOn]);
medTimeUntilContact = nanmedian([flat.wiskContactTimes] - [flat.obsOnTimes]);
medObsOnTime = nanmedian([flat.obsOffTimes] - [flat.obsOnTimes]);
medTrackingTime = nanmedian([flat.lightOnTimes] - [flat.obsOnTimes]);

%% check obstacle tracking...
% use as criterion number of samples where vel deviation is above some
% threshold

% settings
session = '190401_004';
velTime = .01;
obsOnBuffer = .2;
velTolerance = .01;
maxTime = .1;


% initializations
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'obsOnTimes', 'obsOffTimes', 'wheelPositions', 'wheelTimes', 'obsPositions', 'obsTimes', 'targetFs')
wheelVel = getVelocity(wheelPositions, velTime, targetFs);
obsVel = getVelocity(obsPositions, velTime, targetFs);

% get trial data
data = struct();
rowInd = 1;
for i = 1:length(obsOnTimes)
    
    wheelBins = wheelTimes>(obsOnTimes(i)+obsOnBuffer) & wheelTimes<obsOffTimes(i);
    obsBins = obsTimes>obsOnTimes(i) & obsTimes<obsOffTimes(i);
    wheelVelTrial = wheelVel(wheelBins);
    obsVelTrial = obsVel(obsBins);
    times = wheelTimes(wheelBins);
    
    if ~isequal(times, obsTimes(obsBins))
        obsVelTrial = interp1(obsTimes(obsBins), obsVelTrial, times);
    end
    
    data(rowInd).wheelVel = wheelVelTrial;
    data(rowInd).obsVel = obsVelTrial;
    data(rowInd).times = times;
    data(rowInd).goodTracking = nanmean(abs(wheelVelTrial-obsVelTrial));
    rowInd = rowInd + 1;
end

%%


close all; figure; hold on;
plot([data.wheelVel]);
plot([data.obsVel])


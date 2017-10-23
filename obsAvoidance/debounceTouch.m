% function touchDebounced = debounceTouch(touchRaw, touchTimes, breakTimes, obsOnTimes, obsOffTimes)

load('C:\Users\Rick\Google Drive\columbia\obstacleData\sessions\171022_002\run.mat')
load('C:\Users\Rick\Google Drive\columbia\obstacleData\sessions\171022_002\runAnalyzed.mat')
%%

% user settings
minBounce = .01;
touchThresh = 3.5;

% initializations
touchRaw = touch.values;
touchTimes = touch.times;
breakTimes = breaks.times;

% get touch on and off times
onTimes  = touchTimes(logical([0; diff(touchRaw>touchThresh)==1]));
offTimes = touchTimes(logical([0; diff(touchRaw>touchThresh)==-1]));

% ensure first time is on, last is off
onTimes  = onTimes(onTimes<offTimes(end));
offTimes = offTimes(offTimes>onTimes(1));

% for every touch, combine with with subsequent touch if it occurs within minBounce of current touch's off time
for i = 1:(length(onTimes)-1)
    
    if (onTimes(i+1) - offTimes(i)) <= minBounce
        onTimes(i+1) = nan;
        offTimes(i) = nan;
    end
end

sum(isnan(onTimes))
onTimes  = onTimes(~isnan(onTimes));
offTimes = offTimes(~isnan(offTimes));

% for each wheel break, delete subsequent touch on and off signals within the trial



%%
close all; figure;

plot(touchTimes, touchBinary); hold on;
plot(touchTimes, touchDebounced, 'linewidth', 2); hold on;
scatter(breakTimes, ones(size(breakTimes))*.5)

plot(obsTimes, obsPositions / max(obsPositions))

pimpFig
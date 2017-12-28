function [touchDebounced, onTimes, offTimes] = debounceTouch(touchRaw, touchTimes, obsOffTimes, breakTimes)

% takes analog recording of digital touch signal and debounces it into a true digital signal
% returns times at which touches start and end, and continuous touch signal
% adjacent touches are merged together if the interval between their off and on times is less than minBounce
% also, touches ocurring after a wheel break is engaed are removed
%
% inputs        touchRaw:       raw touch signal recorded from spike
%               touchTimes:     sample times for touchRaw
%               obsOffTimes:    times at which the obstacle becomes desengaged (reaches end of track)
%               breakTimes:     times at which the wheel break is engaged
%
% outputs       touchBounced:   continuous debounced touch signal, as described above
%               onTimes:        times at which touches begin
%               offTimes:       times at which touches end


% user settings
minBounce = .02;
touchThresh = 3;
breakTimeThresh = .05; % seconds of touch that activates the wheel break

% get touch on and off times
onTimes  = touchTimes(logical([0; diff(touchRaw>touchThresh)==1]));
offTimes = touchTimes(logical([0; diff(touchRaw>touchThresh)==-1]));

% initialize empty breakTimes vector if not provided
if nargin==3; breakTimes = []; end



if ~isempty(onTimes)

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

    onTimes  = onTimes(~isnan(onTimes));
    offTimes = offTimes(~isnan(offTimes));



    % for each wheel break, delete subsequent touch on and off signals within the trial
    for i = 1:length(breakTimes)
        try
            trialObsOffTime = obsOffTimes( find(obsOffTimes>breakTimes(i), 1, 'first') );
            invalidInds = onTimes>(breakTimes(i)-breakTimeThresh+minBounce) & onTimes<trialObsOffTime;
            onTimes = onTimes(~invalidInds);
            offTimes = offTimes(~invalidInds);
        catch
            keyboard
        end
    end



    % create continuous touch signal
    touchDebounced = zeros(1, length(touchRaw));
    for i = 1:length(onTimes)

        inds = touchTimes>=onTimes(i) & touchTimes<offTimes(i);
        touchDebounced(inds) = 1;

    end
else
    touchDebounced = [];
end




function [sessionVels, isLightOn, obsOnPositions] = getSessionSpeedInfo(session, posInterp, obsPos)

% !!! needs documentation, but generally gets velocities for each trial
% (sessionVels), where each row is a trial // vel is function of position,
% and is interpolated over values specified in posInterp, which is a vector
% of evenly spaced positions // obsPos is the position of the obstacle //
% posInterp is specified relative to obstacle (eg how many meters before
% and after the obs would you like to consider)

% load session data
load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'],...
        'obsPositions', 'obsTimes',...
        'obsOnTimes', 'obsOffTimes',...
        'obsLightOnTimes', 'obsLightOffTimes',...
        'wheelPositions', 'wheelTimes', 'targetFs');
obsPositions = fixObsPositions(obsPositions, obsTimes, obsOnTimes);

% compute velocity
vel = getVelocity(wheelPositions, .5, targetFs);

obsPrePost = [-posInterp(1) posInterp(end)];

% iterate over all trials
sessionVels = nan(length(obsOnTimes), length(posInterp));
isLightOn = false(length(obsOnTimes), 1);
obsOnPositions = nan(length(obsOnTimes), 1);

for j = 1:length(obsOnTimes)

    % locate trial
    obsOnPos = obsPositions( find(obsTimes >= obsOnTimes(j), 1, 'first') );
    obsTime  = obsTimes(find( obsTimes >= obsOnTimes(j) & obsPositions >= obsPos, 1, 'first')); % time at which obstacle reaches obsPos
    obsWheelPos = wheelPositions(find(wheelTimes>=obsTime, 1, 'first')); % position of wheel at moment obstacle reaches obsPos

    % get trial positions and velocities
    trialInds = (wheelPositions > obsWheelPos-obsPrePost(1)) & (wheelPositions < obsWheelPos+obsPrePost(2));
    trialPos = wheelPositions(trialInds);
    trialPos = trialPos - obsWheelPos; % normalize s.t. 0 corresponds to the position at which the obstacle is directly over the wheel
    trialVel = vel(trialInds);

    % remove duplicate positional values
    [trialPos, uniqueInds] = unique(trialPos, 'stable');
    trialVel = trialVel(uniqueInds);

    % interpolate velocities across positional grid and store results
    trialVelInterp = interp1(trialPos, trialVel, posInterp, 'linear');
    sessionVels(j,:) = trialVelInterp;

    % find whether light was on
    isLightOn(j) = min(abs(obsOnTimes(j) - obsLightOnTimes)) < 1; % did the light turn on near whether the obstacle turned on

    % record position at which obstacle turned on
    obsOnPositions(j) = obsPos - obsOnPos;
end
function data = getSpeedAndObsAvoidanceData(sessions)

% creates data structure where each row is a trial for each trial in all of
% sessions // records session, mouse, trialNum, position at which obs
% turned on, trial speed (avg vel while obs is on), whether trial was
% successful, and interpolated speed as a function of position relative to
% obstacle


% settings
obsPrePost = [-.6 .25]; % plot this much before and after the obstacle reaches nose
posRes = .001; % meters // resoultion of positional grid that velocities are computed over
obsContactMinMax = [-.02 .08]; % (meters) only count touches that occur within this range in front of and behind nose


% initializations
posInterp = obsPrePost(1) : posRes : obsPrePost(2); % velocities will be interpolated across this grid of positional values
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')


% collect data
data = struct(); % stores trial data for all sessions
dataInd = 1;

for i = 1:length(sessions)
    
    fprintf('collecting data for session %s...\n', sessions{i});
    
    % get mouse name and experiment name
    sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx']);
    sessionInd = find(strcmp(sessionInfo.session, sessions{i}));

    % load session data
    load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\runAnalyzed.mat'],...
            'obsPositions', 'obsTimes', 'obsPixPositions', 'frameTimeStamps', 'obsOnTimes', 'obsOffTimes',...
            'isLightOn', 'nosePos', 'touchOnTimes', 'touchOffTimes', 'wheelPositions', 'wheelTimes');
    obsPositions = fixObsPositions(obsPositions, obsTimes, obsPixPositions, frameTimeStamps, obsOnTimes, obsOffTimes, nosePos(1));
    wheelVel = getVelocity(wheelPositions, .5, 1/median(diff(wheelTimes)));
    
    % get touch positions and ensure all touches occur when obs is close to mouse
    touchPositions = interp1(obsTimes, obsPositions, touchOnTimes, 'linear');
    validPosBins = touchPositions>obsContactMinMax(1) & touchPositions<obsContactMinMax(2);
    touchOnTimes = touchOnTimes(validPosBins);
    touchOffTimes = touchOffTimes(validPosBins);
    
    % find positions at which obstacles turn on and off
    obsOnPositions = interp1(obsTimes, obsPositions, obsOnTimes, 'linear');
    obsOffPositions = interp1(obsTimes, obsPositions, obsOffTimes, 'linear');
    
    % iterate over all trials
    isObsAvoided = nan(length(obsOnTimes), 1);
    
    for j = 1:length(obsOnTimes)
        
        % find whether and where obstacle was toucheed
        isObsAvoided = ~any(touchOnTimes>obsOnTimes(j) & touchOnTimes<obsOffTimes(j));
        
        % get trial velocity
        if isObsAvoided
            endTime = obsOffTimes(j);
        else
            endTime = touchOnTimes(find(touchOnTimes>obsOnTimes(j),1,'first'));
        end
        wheelInds = find(wheelTimes>=obsOnTimes(j) & wheelTimes<endTime);
        dp = wheelPositions(wheelInds(end)) - wheelPositions(wheelInds(1));
        dt = wheelTimes(wheelInds(end)) - wheelTimes(wheelInds(1));
        avgVel = dp/dt;
        
        
        % get speed as a function of position
        obsWheelPos = interp1(wheelTimes, wheelPositions, obsOnTimes(j)); % wheel position at moment obs turns on
        trialBins = (wheelPositions > obsWheelPos+obsPrePost(1)) & (wheelPositions < obsWheelPos+obsPrePost(2));
        trialPos = wheelPositions(trialBins) - obsWheelPos; % normalize s.t. 0 corresponds to the position at which the obstacle is at the mouse's nose
        trialVel = wheelVel(trialBins); % wheel vel for trial
        
        % remove duplicate positional values
        [trialPos, uniqueInds] = unique(trialPos, 'stable');
        trialVel = trialVel(uniqueInds);

        % interpolate velocities across positional grid and store results
        trialVelInterp = interp1(trialPos, trialVel, posInterp, 'linear');
        
        % store trial info
        data(dataInd).session = sessions{i};
        data(dataInd).mouse = sessionInfo.mouse{sessionInd};
        data(dataInd).experiment = sessionInfo.experiment{sessionInd};
        data(dataInd).trialNum = j;
        data(dataInd).isLightOn = isLightOn(j);
        data(dataInd).isObsAvoided = isObsAvoided;
        data(dataInd).obsOnPositions = obsOnPositions(j);
        data(dataInd).avgVel = avgVel;
        data(dataInd).trialVelInterp = trialVelInterp;
        data(dataInd).trialPosInterp = posInterp;
        dataInd = dataInd+1;
    end
end
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
disp('all done!')
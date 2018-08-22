function data = getSpeedAndObsAvoidanceData(sessions, includeContinuousVelocity)

% creates data structure where each row is a trial for each trial in all of
% sessions // records session, mouse, trialNum, position at which obs
% turned on, trial speed (avg vel while obs is on), whether trial was
% successful, and interpolated speed as a function of position relative to
% obstacle
% if includeContinuousVelocity is true, also computes continuous velocity
% as a function of obs position, which can slow down the function a bit


% settings
obsPrePost = [-.6 .25]; % plot this much before and after the obstacle reaches nose
posRes = .001; % meters // resoultion of positional grid that velocities are computed over


% initializations
posInterp = obsPrePost(1) : posRes : obsPrePost(2); % velocities will be interpolated across this grid of positional values


% collect data for all sessions
data = cell(1,length(sessions));
getDataForSessionHandle = @getDataForSession;
parfor i = 1:length(sessions)
    fprintf('%s: collecting data...\n', sessions{i});
    data{i} = feval(getDataForSessionHandle, sessions{i});
end
data = cat(2,data{:}); % concatenate together data from all sessions 




function sessionData = getDataForSession(session)
    
    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')
    % collect data
    sessionData = struct(); % stores trial data for all sessions
    dataInd = 1;


    fprintf('%s: collecting data...\n', session);

    % get mouse name and experiment name
    sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx']);
    sessionInd = find(strcmp(sessionInfo.session, session));

    % load session data
    load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'],...
            'obsPositions', 'obsTimes', 'obsPixPositions', 'frameTimeStamps', 'obsOnTimes', 'obsOffTimes',...
            'isLightOn', 'nosePos', 'wheelPositions', 'wheelTimes', 'touchesPerPaw');
    load([getenv('OBSDATADIR') 'sessions\' session '\run.mat'], 'breaks');
    obsPositions = fixObsPositions(obsPositions, obsTimes, obsPixPositions, frameTimeStamps, obsOnTimes, obsOffTimes, nosePos(1));
    wheelVel = getVelocity(wheelPositions, .5, 1/median(diff(wheelTimes)));

    % find positions at which obstacles turn on and off
    obsOnPositions = interp1(obsTimes, obsPositions, obsOnTimes, 'linear');

    % iterate over all trials
    for j = 1:length(obsOnTimes)

        % find whether and where obstacle was toucheed
        isWheelBreak = any(breaks.times>obsOnTimes(j) & breaks.times<obsOffTimes(j));

        % get trial velocity
        if isWheelBreak
            endTime = breaks.times(find(breaks.times>obsOnTimes(j),1,'first'));
        else
            endTime = obsOffTimes(j);
        end
        wheelInds = find(wheelTimes>=obsOnTimes(j) & wheelTimes<endTime);
        dp = wheelPositions(wheelInds(end)) - wheelPositions(wheelInds(1));
        dt = wheelTimes(wheelInds(end)) - wheelTimes(wheelInds(1));
        avgVel = dp/dt;


        % get speed as a function of position
        if includeContinuousVelocity
            obsWheelPos = interp1(wheelTimes, wheelPositions, obsOnTimes(j)); % wheel position at moment obs turns on
            trialBins = (wheelPositions > obsWheelPos+obsPrePost(1)) & (wheelPositions < obsWheelPos+obsPrePost(2));
            trialPos = wheelPositions(trialBins) - obsWheelPos; % normalize s.t. 0 corresponds to the position at which the obstacle is at the mouse's nose
            trialVel = wheelVel(trialBins); % wheel vel for trial

            % remove duplicate positional values
            [trialPos, uniqueInds] = unique(trialPos, 'stable');
            trialVel = trialVel(uniqueInds);

            % interpolate velocities across positional grid and store results
            trialVelInterp = interp1(trialPos, trialVel, posInterp, 'linear');
        end

        % get total number of obs touches in trial per paw
        trialBinsTemp = frameTimeStamps>obsOnTimes(j) & frameTimeStamps<obsOffTimes(j);
        trialTouchesPerPaw = touchesPerPaw(trialBinsTemp);
        totalTouchFramesPerPaw = sum(touchesPerPaw(trialBinsTemp,:)>0,1); % get total number of obs touches in trial per paw


        % store trial info
        sessionData(dataInd).session = session;
        sessionData(dataInd).mouse = sessionInfo.mouse{sessionInd};
        sessionData(dataInd).experiment = sessionInfo.experiment{sessionInd};
        sessionData(dataInd).trialNum = j;
        sessionData(dataInd).isLightOn = isLightOn(j);
        sessionData(dataInd).isWheelBreak = isWheelBreak;
        sessionData(dataInd).obsOnPositions = obsOnPositions(j);
        sessionData(dataInd).avgVel = avgVel;
        sessionData(dataInd).trialTouchesPerPaw = trialTouchesPerPaw;
        sessionData(dataInd).totalTouchFramesPerPaw = totalTouchFramesPerPaw;
        if includeContinuousVelocity
            sessionData(dataInd).trialVelInterp = trialVelInterp;
            sessionData(dataInd).trialPosInterp = posInterp;
        end
        dataInd = dataInd+1;
    end
end


warning('on', 'MATLAB:table:ModifiedAndSavedVarnames')
disp('all done!')
end
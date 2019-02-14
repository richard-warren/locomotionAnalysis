function data = getSpeedAndObsAvoidanceData(sessions, sessionInfo, includeContinuousVelocity)

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
velocityDelta = .01;  % compute velocity over this interval

% initializations
posInterp = obsPrePost(1) : posRes : obsPrePost(2); % velocities will be interpolated across this grid of positional values


% collect data for all sessions
data = cell(1,length(sessions));
getDataForSessionHandle = @getDataForSession;
if any(strcmp(sessionInfo.Properties.VariableNames, 'notes'))
    sessionInfo = sessionInfo(:, ~strcmp(sessionInfo.Properties.VariableNames, 'notes'));
end
metaDataFields = sessionInfo.Properties.VariableNames;
metaDataFields = cat(2, metaDataFields, {'sessionNum'});
if ismember('condition', metaDataFields); metaDataFields = cat(2, metaDataFields, {'conditionNum'}); end

for i = 1:length(sessions)
    try
        % get metadata for sessions
        sessionInfoBin = strcmp(sessionInfo.session, sessions{i});
        sessionMetaData = table2struct(sessionInfo(sessionInfoBin,:));
        
        % get sesionNum and conditionNum for mouse
        mouseBins = strcmp(sessionInfo.mouse, sessionMetaData.mouse);
        sessionMetaData.sessionNum = find(strcmp(sessionInfo.session(mouseBins), sessionMetaData.session)); % session num is 1,2,3... for sequential sessions for a given mouse // used to plot performace across days
        if ismember('condition', metaDataFields)
            conditionBins = strcmp(sessionInfo.condition, sessionMetaData.condition);
            sessionMetaData.conditionNum = find(strcmp(sessionInfo.session(mouseBins & conditionBins), sessionMetaData.session)); % first session for condition pre and condition post is 1, second session for condition pre and condition post is 2, etc // used to restrict how many days post lesion to include in analysis
        end
        
        % get session data
        if sessionMetaData.include    
            data{i} = feval(getDataForSessionHandle, sessions{i}, sessionMetaData);
        else
            fprintf('%s: skipped\n', sessions{i})
        end
        
    catch
        fprintf('%s: unable to analyze session!\n', sessions{i});
    end
end
data = cat(2,data{:}); % concatenate together data from all sessions 




function sessionData = getDataForSession(session, sessionMetaData)
    
    % collect data
    fprintf('%s: collecting data...\n', session);
    sessionData = struct(); % stores trial data for all sessions
    dataInd = 1;

    % load session data
    load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'],...
            'obsPositions', 'obsTimes', 'obsPixPositions', 'frameTimeStamps', 'obsOnTimes', 'obsOffTimes',...
            'isLightOn', 'nosePos', 'wheelPositions', 'wheelTimes', 'wiskContactPositions', ...
            'touches', 'touchesPerPaw', 'touchClassNames', 'bodyAngles', 'obsPositionsFixed');
    load([getenv('OBSDATADIR') 'sessions\' session '\run.mat'], 'breaks');
    wheelVel = getVelocity(wheelPositions, velocityDelta, 1/median(diff(wheelTimes)));

    % find positions at which obstacles turn on and off
    obsOnPositions = interp1(obsTimes, obsPositionsFixed, obsOnTimes, 'linear');

    % iterate over all trials
    for j = 1:length(obsOnTimes)

        % find whether and where obstacle was toucheed
        isWheelBreak = any(breaks.times>obsOnTimes(j) & breaks.times<obsOffTimes(j));

        % get trial velocity and body angle
        if isWheelBreak
            endTime = breaks.times(find(breaks.times>obsOnTimes(j),1,'first'));
        else
            endTime = obsOffTimes(j);
        end
        wheelInds = find(wheelTimes>=obsOnTimes(j) & wheelTimes<endTime);
        dp = wheelPositions(wheelInds(end)) - wheelPositions(wheelInds(1));
        dt = wheelTimes(wheelInds(end)) - wheelTimes(wheelInds(1));
        avgVel = dp/dt;
        avgAngle = nanmedian(bodyAngles(frameTimeStamps>obsOnTimes(j) & frameTimeStamps<endTime));


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
        trialTouches = touches(trialBinsTemp);
        trialTouchesPerPaw = touchesPerPaw(trialBinsTemp,:);
        totalTouchFramesPerPaw = sum(touchesPerPaw(trialBinsTemp,:)>0,1); % get total number of obs touches in trial per paw


        % session meta data
        % session metadata (determined from sessionInfo table)
        for k = metaDataFields
            sessionData(dataInd).(k{1}) = sessionMetaData.(k{1});
        end

        % trial metadata
        sessionData(dataInd).trialNum = j;
        sessionData(dataInd).isLightOn = isLightOn(j);
        sessionData(dataInd).isWheelBreak = isWheelBreak;
        sessionData(dataInd).obsOnPositions = obsOnPositions(j);
        sessionData(dataInd).avgVel = avgVel;
        sessionData(dataInd).avgAngle = avgAngle;
        sessionData(dataInd).trialTouches = trialTouches;
        sessionData(dataInd).trialTouchesPerPaw = trialTouchesPerPaw;
        sessionData(dataInd).totalTouchFramesPerPaw = totalTouchFramesPerPaw;
        sessionData(dataInd).touchClassNames = touchClassNames;
        if includeContinuousVelocity
            sessionData(dataInd).trialVelInterp = trialVelInterp;
            sessionData(dataInd).trialPosInterp = posInterp;
            sessionData(dataInd).wiskContactPositions = wiskContactPositions;
        end
        dataInd = dataInd+1;
    end
end

disp('all done!')
end
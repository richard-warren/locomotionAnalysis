function expData = getExperimentData(sessionInfo, vars, oldData)

% creates nested struct containing data for all mice, sessions, trials, and
% paws // at each level in this heirarchy vars can be computed // add
% sanity checks for all dvs (eg step over hgt cant be less than hgt of asdf
% obstacle...) // dft fcns for getting session, mouse, trial vars, etc...

% todo: store all sessionInfo columns as session variables automatically //
% change conditionNum and sessionNum to colums in spreadsheet to handle day
% discontinuities, or you must deal with dates...

% 'g' stores global variables, allowing things to be passed around!
% some of this code assumes the first six characters of the session name
% encode the date in the format: YYMMDD


% settings
metadata = {'touchThresh', 'speedTime', 'preObsLim', 'clearanceBuffer', 'velVsPositionX', 'velContinuousAtContactX'};  % these parameters will be stored as experiment metadata
g.touchThresh = 5;  % successful trials have fewer than touchThresh frames where paw is in contact with obs
g.speedTime = .01;  % (s) compute velocity over this interval
g.preObsLim = .008;  % (m) compute paw height this many meters before it reaches obstacle x postion
g.clearanceBuffer = .000;  % (m) trials are excluded in which paw height is less that obsHeight - pawClearnceBuffer at the moment it reaches the x position of the obstacle

g.velVsPositionPrePost = [-1.5 .4]; % (m) positions before and after whisker contact to collect velocity data
g.velVsPositionRes = 500; % (tics) how many samples in the x grid

g.velContinuousAtContactPrePost = [-1 1]; % (s) how many seconds before and after wisk contact to compute velocity
g.velContinuousAtContactRes = 500; % (tics) how many samples in the x grid


% initialiations
if ischar(sessionInfo) % if sessionInfo is a string, then it contains the name of a single session, and we will get sesionInfo for this session automatically
    sessionName = sessionInfo;
    sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'sessionInfo.xlsx'), 'Sheet', 'sessions');
    sessionInfo = sessionInfo(strcmp(sessionInfo.session, sessionName),:);
end
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);  % remove not included sessions and empty lines
g.sessionInfo = sessionInfo;
g.velVsPositionX = linspace(g.velVsPositionPrePost(1), g.velVsPositionPrePost(2), g.velVsPositionRes);
g.velContinuousAtContactX = linspace(g.velContinuousAtContactPrePost(1), g.velContinuousAtContactPrePost(2), g.velContinuousAtContactRes);


mouseVars = {};
sessionVars = {'condition', 'side', 'brainRegion', 'conditionNum', 'sessionNum', 'whiskers'};
trialVars = {'obsOnTimes', 'obsOffTimes', 'obsOnPositions', 'obsOffPositions', 'velContinuousAtContact', 'velVsPosition', 'isLightOn', ...
             'isWheelBreak', 'obsHgt', 'isTrialSuccess', 'trialVel', 'velAtWiskContact', 'firstModPaw', ...
             'trialAngle', 'trialAngleContra', 'angleAtWiskContact', 'angleAtWiskContactContra', ...
             'wiskContactPosition', 'wiskContactTimes', 'lightOnTimes', 'isContraFirst', 'isBigStep', 'isModPawContra', ...
             'tailHgt', 'tailHgtAtWiskContact', 'modPawDistanceToObs', 'modPawPredictedDistanceToObs', 'velContinuousAtContact', ...
             'modPawKin', 'modPawKinInterp', 'preModPawKin', 'preModPawKinInterp', 'modPawDeltaLength', 'preModPawDeltaLength', ...
             'sensoryCondition', 'contactInd', 'contactIndInterp', 'trialDuration', 'touchFrames', 'modPawOnlySwing'};
pawVars = {'isContra', 'isFore', 'isLeading', 'isPawSuccess', 'stepOverMaxHgt', 'preObsHgt', 'controlPreObsHgt', 'controlStepHgt', 'noObsStepHgt', ...
           'stepOverStartingDistance', 'stepOverEndingDistance', 'stepOverKinInterp', 'controlStepKinInterp', ...
           'isValidZ', 'preObsKin', 'xDistanceAtPeak', 'stepOverLength', 'preStepOverLength', 'prePreStepOverLength', 'controlStepLength', ...
           'isVentralContact', 'isDorsalContact', 'numTouchFrames', 'stepType', 'distanceToObs', 'anyTouchFrames'};

% compute only requested vars
if isequal(vars, 'all'); vars = cat(2, sessionVars, trialVars, pawVars); end
mouseVars = mouseVars(ismember(mouseVars, vars));
sessionVars = sessionVars(ismember(sessionVars, vars));
trialVars = trialVars(ismember(trialVars, vars));
pawVars = pawVars(ismember(pawVars, vars));



% get mouse data
g.mice = unique(g.sessionInfo.mouse);
g.expData = struct('mouse', g.mice);
for mouseVar = 1:length(mouseVars)
    temp = getVar(mouseVars{mouseVar});
    [g.expData(1:length(g.mice)).(mouseVars{mouseVar})] = temp{:};
end

% loop over mice
disp('getting experiment data...')
for mouse = 1:length(g.mice)
    g.mouse = mouse;
    
    % get session data
    g.sessions = g.sessionInfo.session(strcmp(g.sessionInfo.mouse, g.mice{mouse}));
    g.expData(mouse).sessions = struct('session', g.sessions);  % create nested 'session' struct for mouse
    for sessionVar = 1:length(sessionVars)
        try
        temp = getVar(sessionVars{sessionVar});
        [g.expData(mouse).sessions(1:length(g.sessions)).(sessionVars{sessionVar})] = temp{:};
        catch; end
    end
    
    % loop over sessions
    for session = 1:length(g.sessions)
        g.session = session;
        
        % check if session already exists in oldData
        sesBin = false;
        if exist('oldData', 'var')
            mouseBin = strcmp({oldData.data.mouse}, g.mice{mouse});
            if any(mouseBin)
                sesBin = strcmp({oldData.data(mouseBin).sessions.session}, g.sessions{session});
            end
        end
        
        % if session exists in old data, copy it over
        if any(sesBin)
            fprintf('%s: copying trials from previous data...\n', g.sessions{session})
            g.expData(mouse).sessions(session).trials = oldData.data(mouseBin).sessions(sesBin).trials;
        
        % otherwise, compute de novo
        else
            % load data from session, or compute if necessary
            if exist(fullfile(getenv('OBSDATADIR'), 'sessions', g.sessions{session}, 'kinData.mat'), 'file')
                g.sesKinData = load(fullfile(getenv('OBSDATADIR'), 'sessions', g.sessions{session}, 'kinData.mat'), 'kinData');
            else
                try
                    g.sesKinData = getKinematicData(g.sessions{session});
                catch
                    fprintf('%s: problem with getKinematicData!', g.sessions{session});
                end
            end
            g.sesKinData = g.sesKinData.kinData;
            g.sesKinInds = find([g.sesKinData.isTrialAnalyzed]);
            g.sesData = load(fullfile(getenv('OBSDATADIR'), 'sessions', g.sessions{session}, 'runAnalyzed.mat'));
            warning('off', 'MATLAB:load:variableNotFound');
            g.sesSpikeData = load(fullfile(getenv('OBSDATADIR'), 'sessions', g.sessions{session}, 'run.mat'), 'breaks', 'stimulus');
            g.wheelVel = getVelocity(g.sesData.wheelPositions, g.speedTime, g.sesData.targetFs);
            g.isValidZ = getSessionValidZ(g.sesKinData, g.sesData.obsHeights/g.sesData.targetFs);

            % get size of kin data entries
            g.locationsSmps = size(g.sesKinData(g.sesKinInds(1)).modifiedLocations{1}, 3);
            g.locationsInterpSmps = size(g.sesKinData(g.sesKinInds(1)).modifiedLocationsInterp{1}, 3);
            g.numControlSteps = size(g.sesKinData(g.sesKinInds(1)).controlLocations{1},1);

            % get trial data
            g.expData(mouse).sessions(session).trials = struct('trial', num2cell(1:length(g.sesKinData)));  % create nested 'trial' struct for mouse
            for trialVar = 1:length(trialVars)
                try
                temp = getVar(trialVars{trialVar});
                [g.expData(mouse).sessions(session).trials(1:length(g.sesKinData)).(trialVars{trialVar})] = temp{:};
                catch; end
            end

            % loop over trials
            isStepOverAnalyzed = false(length(g.sesKinData),4);
            for trial = 1:length(g.sesKinData)
                g.trial = trial;

                % get paw data
                g.expData(mouse).sessions(session).trials(trial).paws = struct('paw', {1,2,3,4}, 'pawName', {'LH', 'LF', 'RF', 'RH'});  % create nested 'paw' struct for mouse
                for pawVar = 1:length(pawVars)
                    try
                    temp = getVar(pawVars{pawVar});
                    [g.expData(mouse).sessions(session).trials(trial).paws(1:4).(pawVars{pawVar})] = temp{:};
                    catch; end
                end
                
                % record whether each step is successfully analyzed
                if ismember('stepOverKinInterp', vars)
                    isStepOverAnalyzed(trial,:) = cellfun(@(x) ~all(isnan(x(:))), {g.expData(mouse).sessions(session).trials(trial).paws.stepOverKinInterp});
                end
            end

            stepOverAnalyzedRates = num2cell(mean(isStepOverAnalyzed,1));
            if ~isempty(pawVars)
                fprintf('%s: analysis success rates -> paw1: %.2f, paw2: %.2f, paw3: %.2f, paw4: %.2f \n', ...
                    g.sessions{session}, stepOverAnalyzedRates{:})
            end
        end
    end
end

% save experiment metadata
expData = g.expData;
dataTemp = expData; clear expData;  % nest expData within itself
expData.data = dataTemp; clear dataTemp;
for m = 1:length(metadata); expData.(metadata{m}) = g.(metadata{m}); end % add one field per metadatum










% ---------
% FUNCTIONS
% ---------

function var = getVar(dvName) % sessionInfo, expData, mice, mouse, sessions, session, trial, sesKinData, sesData, wheelVel, isValidZ, breakTimes, touchThresh, sesKinInds, locationsInterpSmps, preObsLim
    
    switch dvName
        
        % session variables
        % -----------------
        case 'condition'  % read from sessionInfo column
            var = getTableColumn('condition');
            
        case 'side'  % read from sessionInfo column
            var = getTableColumn('side');
            
        case 'brainRegion'  % read from sessionInfo column
            var = getTableColumn('brainRegion');
            
        case 'conditionNum'  % which day of a particular condition is the session, e.g. conditionNum=3 means this is the third session for this condition
            % note: excluded sessions are ignored in conditionNum
            % computation! also, assumes sessions are in temporal order
            % (oldest to most recent) in sessionInfo
            sessionInfoSub = g.sessionInfo(strcmp(g.sessionInfo.mouse, g.expData(g.mouse).mouse),:);
            var = num2cell(nan(1,height(sessionInfoSub)));
            conditions = unique(sessionInfoSub.condition);
            for i = conditions'
                bins = strcmp(sessionInfoSub.condition, i{1});
                var(bins) = num2cell(1:sum(bins));
            end
            var = var(logical(sessionInfoSub.include));
            
        case 'sessionNum'  % day number, e.g. sessionNum=3 means this is the third session recorded in sessionInfo
            % note: excluded sessions are signored in sessionNum computation!
            allMouseSessions = g.sessionInfo.session(strcmp(g.sessionInfo.mouse, g.mice{g.mouse})); % all sessions, including those where .include==false
            [~, inds] = intersect(allMouseSessions, g.sessions, 'stable');
            var = num2cell(inds);
            
        case 'whiskers'
            var = getTableColumn('whiskers');
        
            
            

        % trial variables
        % ---------------
        
        case 'obsOnTimes'  % times at which obstacle turns on
            var = num2cell([g.sesData.obsOnTimes]);
        
        case 'obsOffTimes'  % times at which obstacle turns on
            var = num2cell([g.sesData.obsOffTimes]);
        
        case 'obsOnPositions'  % position of the obstacle relative to mouse nose at the moment it turns on
            var = num2cell(interp1(g.sesData.obsTimes, g.sesData.obsPositionsFixed, g.sesData.obsOnTimes, 'linear'));
            
        case 'obsOffPositions'  % position of the obstacle relative to mouse nose at the momentit turns off
            var = num2cell(interp1(g.sesData.obsTimes, g.sesData.obsPositionsFixed, g.sesData.obsOffTimes, 'linear'));
        
        case 'velContinuousAtContact'  % continuous velocity vector surrouding moment of whisker contact
            times = g.velContinuousAtContactX;
            var = repmat({nan(1, length(times))}, 1, length(g.sesKinData));
            
            for i = g.sesKinInds
                bins = g.sesData.wheelTimes>(g.sesData.obsOnTimes(i)-range(times)) & ...
                    g.sesData.wheelTimes<(g.sesData.obsOffTimes(i)+range(times));  % add range(times) as a buffer in case the desired times fall outside the range of the obstacle being on
                var{i} = interp1(g.sesData.wheelTimes(bins), g.wheelVel(bins), times+g.sesData.wiskContactTimes(i));
            end
            
        case 'velVsPosition'  % continuous velocity as a function of position of obstacle relative to nose
            var = repmat({nan(1, length(g.velVsPositionX))}, 1, length(g.sesKinData));
            
            for i = g.sesKinInds
                obsAtNoseTime = g.sesData.obsTimes(find(g.sesData.obsPositionsFixed>=0 & ...
                                g.sesData.obsTimes>g.sesData.obsOnTimes(i) & ...
                                g.sesData.obsTimes<g.sesData.obsOffTimes(i), 1, 'first'));
                if ~isempty(obsAtNoseTime)
                    obsAtNosePos = g.sesData.wheelPositions(find(g.sesData.wheelTimes>obsAtNoseTime,1,'first'));  % position of the WHEEL when obstacle reaches nose
                    trialBins = (g.sesData.wheelPositions > (obsAtNosePos+g.velVsPositionPrePost(1))) & ...
                                (g.sesData.wheelPositions < (obsAtNosePos+g.velVsPositionPrePost(2)));
                    trialPos = g.sesData.wheelPositions(trialBins) - obsAtNosePos;  % normalize s.t. 0 corresponds to the position at which the obstacle is at the mouse's nose
                    trialVel = g.wheelVel(trialBins);  % wheel vel for trial

                    % remove duplicate positional values (would be better to average all values in a particular bin)
                    [trialPos, uniqueInds] = unique(trialPos, 'stable');
                    trialVel = trialVel(uniqueInds);

                    % interpolate velocities across positional grid and store results
                    var{i} = interp1(trialPos, trialVel, g.velVsPositionX, 'linear');
                end
            end
            
        case 'isLightOn'  % whether obstacle light is on for each trial
            var = num2cell(g.sesData.isLightOn);
            
        case 'isWheelBreak'  % whether wheel break is engaged for each trial
            var = num2cell(g.sesData.isWheelBreak);
            
        case 'obsHgt'  % height of obstacle, in meters
            var = num2cell(g.sesData.obsHeights/1000);  % convert back to meters
        
        case 'isTrialSuccess'  % whether paw is touching obstacle in fewer than 'touchThresh' frames for each trial
            var = cell(1,length(g.sesKinData));
            for i = 1:length(g.sesKinData)
                var{i} = sum(any(g.sesData.touchesPerPaw(g.sesKinData(i).trialInds,:),2)) < g.touchThresh;
            end
            
        case 'trialVel'  % average trial velocity, excluding periods after wheel breaks
            var = avgSignalPerTrial(g.wheelVel, g.sesData.wheelTimes);
        
        case 'velAtWiskContact'  % running speed at moment of whisker contact
            var = num2cell(interp1(g.sesData.wheelTimes, g.wheelVel, g.sesData.wiskContactTimes));
            
        case 'firstModPaw'  % first modified paw identity
            var = num2cell(nan(1,length(g.sesKinData)));
            var([g.sesKinData.isTrialAnalyzed]) = num2cell([g.sesKinData.firstModPaw]);
            
        case 'trialAngle'  % average body angle for trial, excluding periods after wheel breaks
            var = avgSignalPerTrial(g.sesData.bodyAngles, g.sesData.frameTimeStamps);
            
        case 'trialAngleContra'  % average body angle for trial wrt 'side', excluding periods after wheel breaks
            side = g.expData(g.mouse).sessions(g.session).side;
            if strcmp(side, 'left')
                var = num2cell(-[g.expData(g.mouse).sessions(g.session).trials.trialAngle]);
            elseif strcmp(side, 'right')
                var = num2cell([g.expData(g.mouse).sessions(g.session).trials.trialAngle]);
            end
        
        case 'angleAtWiskContact'  % body angle at moment of whisker contact
            bins = ~isnan(g.sesData.frameTimeStamps);
            var = num2cell(interp1(g.sesData.frameTimeStamps(bins), g.sesData.bodyAngles(bins), g.sesData.wiskContactTimes));
            
        case 'angleAtWiskContactContra'  % body angle at moment of whisker contact wrt 'side'
            side = g.expData(g.mouse).sessions(g.session).side;
            if strcmp(side, 'left')
                var = num2cell(-[g.expData(g.mouse).sessions(g.session).trials.angleAtWiskContact]);
            elseif strcmp(side, 'right')
                var = num2cell([g.expData(g.mouse).sessions(g.session).trials.angleAtWiskContact]);
            end
            
        case 'wiskContactPosition'  % position of obstacle relative to nose at moment of contact (more negative numbers mean further away from the nose at whisker contact)
            var = num2cell(nan(1,length(g.sesKinData)));
            var([g.sesKinData.isTrialAnalyzed]) = num2cell([g.sesKinData.wiskContactPositions]);
            
        case 'wiskContactTimes'  % times of whisker contact
            var = num2cell(nan(1,length(g.sesKinData)));
            var([g.sesKinData.isTrialAnalyzed]) = num2cell([g.sesKinData.wiskContactTimes]);
            
        case 'lightOnTimes'  % times at which obstacle light turns on
            var = {g.sesData.obsLightOnTimes};
            
        case 'isContraFirst'  % whether the paw contralateral to 'side' is the first paw over the obstacle
            side = g.expData(g.mouse).sessions(g.session).side;
            sequences = cat(2, g.sesKinData.pawOverSequence);
            if strcmp(side, 'left')
                var = num2cell(ismember(sequences(1,:), [3 4]));
            elseif strcmp(side, 'right')
                var = num2cell(ismember(sequences(1,:), [1 2]));
            end
            
        case 'isBigStep'  % whether this is a 'big step' or 'little step' trial
            var = num2cell(nan(1,length(g.sesKinData)));
            for i = g.sesKinInds
                var{i} = size(g.sesKinData(i).modifiedLocations{g.sesKinData(i).firstModPaw},1)==1;
            end
            
        case 'isModPawContra'  % whether the 'first modified paw' (the forepaw in the air at the moment of whisker contact) is contralateral to 'side'
            side = g.expData(g.mouse).sessions(g.session).side;
            if strcmp(side, 'left') || strcmp(side, 'right')
                var = num2cell(false(1,length(g.sesKinData)));
                if strcmp(side, 'right'); contraPaw = 2; else; contraPaw = 3; end
                var([g.sesKinData.isTrialAnalyzed]) = num2cell([g.sesKinData.firstModPaw]==contraPaw);
            end
            
        case 'tailHgt'  % average tail height for trial, excluding periods after wheel breaks
            tailHgts = nan(size(g.sesData.frameTimeStamps));
            for i = g.sesKinInds
                tailHgts(g.sesKinData(i).trialInds) = g.sesKinData(i).locationsTail(:,3,1);
            end
            var = avgSignalPerTrial(tailHgts, g.sesData.frameTimeStamps);
            
        case 'tailHgtAtWiskContact'  % height of base of tail at moment of whisker contact
            var = num2cell(nan(1,length(g.sesKinData)));
            for i = g.sesKinInds
                t = g.sesData.frameTimeStamps(g.sesKinData(i).trialInds);
                z = g.sesKinData(i).locationsTail(:,3,1);
                bins = ~isnan(t) & ~isnan(z);
                var{i} = interp1(t(bins), z(bins), g.sesKinData(i).wiskContactTimes);
            end
            
        case 'modPawDistanceToObs'  % distance of mod paw to obstacle at END of first mod step
            var = num2cell(false(1,length(g.sesKinData)));
            for i = g.sesKinInds
                var{i} = g.sesKinData(i).modifiedLocations{g.sesKinData(i).firstModPaw}(1,1,end);
            end
            
        case 'modPawPredictedDistanceToObs'  % predicted distance of mod paw to obstacle at END of first mod step (where would the paw have landed relative to the obstacle if there were no behavioral modifications?)
            var = num2cell(false(1,length(g.sesKinData)));
            for i = g.sesKinInds
                xStart = g.sesKinData(i).modifiedLocations{g.sesKinData(i).firstModPaw}(1,1,1); % first x val of first mod step for mod paw
                predictedLength = g.sesKinData(i).modPredictedLengths(1,g.sesKinData(i).firstModPaw);
                var{i} = xStart + predictedLength;
            end
            
        case 'modPawKin'  % kinematics of first modified step for first modified paw
            var = repmat({nan(3,g.locationsSmps)},1,length(g.sesKinData));
            isBigStep = [g.expData(mouse).sessions(session).trials.isBigStep];
            for i = g.sesKinInds
                var{i} = squeeze(g.sesKinData(i).modifiedLocations{g.sesKinData(i).firstModPaw}(1,:,:));
                if ~g.isValidZ(i,g.sesKinData(i).firstModPaw) && isBigStep(i)  % if step doesn't pass over obstacle but was supposed to pass over obstacle (i.e. is big step)
                    var{i} = nan(3,g.locationsSmps);
                end
            end
            
        case 'modPawKinInterp'  % interpolated kinematics of first modified step for first modified paw
            var = repmat({nan(3,g.locationsInterpSmps)},1,length(g.sesKinData));
            isBigStep = [g.expData(mouse).sessions(session).trials.isBigStep];
            for i = g.sesKinInds
                var{i} = squeeze(g.sesKinData(i).modifiedLocationsInterp{g.sesKinData(i).firstModPaw}(1,:,:));
                if ~g.isValidZ(i,g.sesKinData(i).firstModPaw) && isBigStep(i)  % if step doesn't pass over obstacle but was supposed to pass over obstacle (i.e. is big step)
                    var{i} = nan(3,g.locationsInterpSmps);
                end
            end
            
        case 'preModPawKin'
            % kinematics of steps (plural!) preceding first modified step for first modified paw
            % note: this matrix contains MULTIPLE steps, and has 3 dimensions as a results ([step X dimension (xyz) X time])
            var = repmat({nan(g.numControlSteps,3,g.locationsSmps)}, 1, length(g.sesKinData));
            for i = g.sesKinInds
                var{i} = squeeze(g.sesKinData(i).controlLocations{g.sesKinData(i).firstModPaw});
            end
            
        case 'preModPawKinInterp'  % interpolated kinematics of step preceding first modified step for first modified paw
            var = repmat({nan(3,g.locationsInterpSmps)},1,length(g.sesKinData));
            for i = g.sesKinInds
                var{i} = squeeze(g.sesKinData(i).controlLocationsInterp{g.sesKinData(i).firstModPaw}(end,:,:));
            end
            
        case 'modPawDeltaLength'  % change in length of first modified step of first modified paw (relative to what we predict would have happened if there were no behavioral modifications)
            var = num2cell(false(1,length(g.sesKinData)));
            for i = g.sesKinInds
                predictedLength = g.sesKinData(i).modPredictedLengths(1,g.sesKinData(i).firstModPaw);
                actualLength = g.sesKinData(i).modifiedSwingLengths(1,g.sesKinData(i).firstModPaw);
                var{i} = actualLength - predictedLength;
            end
            
        case 'preModPawDeltaLength'  % change in length of step preceding first modified step of first modified paw (serves as control for modPawDeltaLength)
            var = num2cell(false(1,length(g.sesKinData)));
            for i = g.sesKinInds
                predictedLength = g.sesKinData(i).controlPredictedLengths(end,g.sesKinData(i).firstModPaw);
                actualLength = g.sesKinData(i).controlSwingLengths(end,g.sesKinData(i).firstModPaw);
                var{i} = actualLength - predictedLength;
            end
            
        case 'sensoryCondition'  % encodes whether animal has light + whiskers ('LW'), whiskers only ('W'), light only ('L'), or neither ('-')
            if ismember('whiskers', g.sessionInfo.Properties.VariableNames)
                hasWhiskers = ~strcmp(g.sessionInfo.whiskers{strcmp(g.sessionInfo.session, g.sessions{g.session})}, 'none');
            else
                hasWhiskers = true;
            end
            if hasWhiskers; conditions = {'W', 'WL'}; else; conditions = {'-', 'L'}; end
            var = conditions([g.sesData.isLightOn]+1);
            
        case 'contactInd'  % ind in modPawKin at which whiskers contact obstacle
            var = num2cell(nan(1,length(g.sesKinData)));
            for i = g.sesKinInds
                var{i} = g.sesKinData(i).pawObsPosInd(g.sesKinData(i).firstModPaw); 
            end
            
        case 'contactIndInterp'  % ind in modPawKinInterp at which whiskers contact obstacle
            var = num2cell(nan(1,length(g.sesKinData)));
            for i = g.sesKinInds
                var{i} = g.sesKinData(i).pawObsPosIndInterp(g.sesKinData(i).firstModPaw); 
            end
            
        case 'trialDuration'  % duration of trial (from obstalce on to obstacle off time)
            var = num2cell(g.sesData.obsOffTimes - g.sesData.obsOnTimes);
            
        case 'touchFrames'  % number of frames in trial where paw is contacting obstacle
            var = cell(1,length(g.sesKinData));
            for i = 1:length(g.sesKinData)
                var{i} = sum(any(g.sesData.touchesPerPaw(g.sesKinData(i).trialInds,:),2));
            end
            
        case 'modPawOnlySwing'  % whether the mod paw is the only one in swing at the moment of contact (see getKinematicData for formal definition of firstModPaw)
            var = nan(size(g.sesKinData));
            var([g.sesKinData.isTrialAnalyzed]) = ...
                xor([g.sesKinData.isRightSwingAtContact], [g.sesKinData.isLeftSwingAtContact]);
            var = num2cell(var);
            
            
            

        % paw variables
        % -------------
        case 'isContra'  % whether paw is contralateral to 'side'
            side = g.expData(g.mouse).sessions(g.session).side;
            if strcmp(side, 'left')
                var = num2cell(logical([0 0 1 1]));
            elseif strcmp(side, 'right')
                var = num2cell(logical([1 1 0 0]));
            end
            
        case 'isFore'  % whether paw is a forelimb
            var = num2cell(logical([0 1 1 0]));
        
        case 'isLeading'  % whether paw is leading (first to get over obstacle relative to other fore/hind limb)
            var = num2cell(logical([1 1 0 0]));  % start assuming left is leading
            seq = g.sesKinData(g.trial).pawOverSequence;
            if find(seq==4)<find(seq==1); var([1,4]) = var([4,1]); end  % if right hind crosses before left hind
            if find(seq==3)<find(seq==2); var([2,3]) = var([3,2]); end  % if right fore crosses before left fore
            
        case 'isPawSuccess'  % whether paw is in contact with obstacle for fewer than touchThresh frames in trial
            var = num2cell(sum(g.sesData.touchesPerPaw(g.sesKinData(g.trial).trialInds,:)~=0,1) < g.touchThresh);
        
        case 'stepOverMaxHgt'  % maximum height of paw on step over obstacle
            var = num2cell(nan(1,4));
            if g.sesKinData(g.trial).isTrialAnalyzed
                var(g.isValidZ(g.trial,:)) = num2cell(cellfun(@(x) max(x(end,3,:)), ...
                    g.sesKinData(g.trial).modifiedLocationsInterp(g.isValidZ(g.trial,:))));
                if any(cellfun(@(x) x<0, var(g.isValidZ(g.trial,:)))); disp('trial included with NEGATIVE max height... what the fuck?!'); end
            end
            
        case 'preObsHgt'  % height of paw when the step over is preObsLim in front of osbtacle (use to gauge height shaping at a moment when the paw has not had the opportunity to contact the osbatcle)
            var = num2cell(nan(1,4));
            if g.sesKinData(g.trial).isTrialAnalyzed
                for i = find(g.isValidZ(g.trial,:))
                    ind = find(g.sesKinData(g.trial).modifiedLocationsInterp{i}(end,1,:) > (-g.preObsLim), 1, 'first');
                    if ind>1; var{i} = g.sesKinData(g.trial).modifiedLocationsInterp{i}(end,3,ind); end
                end
            end
            
        case 'controlPreObsHgt'  % for final control step, find height at same index in step over at which paw reaches preObsLim in front of obstalce // this is intended to allow comparison of preObsHgt to a similar position in control steps
            var = num2cell(nan(1,4));
            if g.sesKinData(g.trial).isTrialAnalyzed
                for i = find(g.isValidZ(g.trial,:))
                    ind = find(g.sesKinData(g.trial).modifiedLocationsInterp{i}(end,1,:) > (-g.preObsLim), 1, 'first');
                    if ind>1; var{i} = g.sesKinData(g.trial).controlLocationsInterp{i}(end,3,ind); end
                end
            end
            
        case 'noObsStepHgt'  % avg step hgt for steps before the obstalce turns on
            var = num2cell(nan(1,4));
            if g.sesKinData(g.trial).isTrialAnalyzed
                var = cellfun(@(x) nanmean(max(squeeze(x(:,3,:)),[],2)), ...  % avg the max z for each noObsStep
                    g.sesKinData(g.trial).noObsLocationsInterp, 'UniformOutput', false);
            end
            
        case 'controlStepHgt'  % avg step hgt for steps before first modified step
            var = num2cell(nan(1,4));
            if g.sesKinData(g.trial).isTrialAnalyzed
                var = cellfun(@(x) nanmean(max(squeeze(x(:,3,:)),[],2)), ... % avg the max z for each noObsStep
                    g.sesKinData(g.trial).controlLocationsInterp, 'UniformOutput', false);
            end
            
        case 'stepOverStartingDistance'  % distance of paw to obs at start of step over obs
            var = num2cell(nan(1,4));
            if g.sesKinData(g.trial).isTrialAnalyzed
                var = num2cell(cellfun(@(x) x(end,1,1), g.sesKinData(g.trial).modifiedLocations));
            end
            
        case 'stepOverEndingDistance'  % distance of paw to obs at end of step over obs
            var = num2cell(nan(1,4));
            if g.sesKinData(g.trial).isTrialAnalyzed
                var = num2cell(cellfun(@(x) x(end,1,end), g.sesKinData(g.trial).modifiedLocationsInterp));
            end
        
        case 'stepOverKinInterp'  % interpolated kinematics of step over obstacle
            var = repmat({nan(3,g.locationsInterpSmps)},1,4);
            if g.sesKinData(g.trial).isTrialAnalyzed
                var(g.isValidZ(g.trial,:)) = cellfun(@(x) squeeze(x(end,:,:)), ...
                    g.sesKinData(g.trial).modifiedLocationsInterp(g.isValidZ(g.trial,:)), 'UniformOutput', false);
            end
            
        case 'controlStepKinInterp' % interpolated kinematics of last control step
            var = repmat({nan(3,g.locationsInterpSmps)},1,4);
            if g.sesKinData(g.trial).isTrialAnalyzed
                var(g.isValidZ(g.trial,:)) = cellfun(@(x) squeeze(x(end,:,:)), ...
                    g.sesKinData(g.trial).controlLocationsInterp(g.isValidZ(g.trial,:)), 'UniformOutput', false);
            end
            
        case 'isValidZ'
            % checks whether step over obstacle passes under obstacle, indicating a tracking error
            % other variables that rely on z coordinates should check whether isValidZ
            var = num2cell(g.isValidZ(g.trial,:));
            
        case 'preObsKin'  % interpolated kinematics of step over obs, but truncated at the moment paw is preObsLim in from the obstacle
            var = repmat({nan(3,g.locationsInterpSmps)},1,4);
            if g.sesKinData(g.trial).isTrialAnalyzed
                for i = find(g.isValidZ(g.trial,:))
                    ind = find(g.sesKinData(g.trial).modifiedLocations{i}(end,1,:) >- (g.preObsLim), 1, 'first');
                    kin = squeeze(g.sesKinData(g.trial).modifiedLocations{i}(end,:,1:ind));
                    if ind>1
                        for j = 1:3 % x, y, z
                            var{i}(j,:) = interp1(1:size(kin,2), kin(j,:), linspace(1,size(kin,2), g.locationsInterpSmps), 'linear');
                        end 
                    end
                end
            end
            
        case 'xDistanceAtPeak'  % horizontal (AP) distance of paw to obstacle when it is at max height
            var = num2cell(nan(1,4));
            if g.sesKinData(g.trial).isTrialAnalyzed
                for i = find(g.isValidZ(g.trial,:))
                    [~, maxInd] = max(g.sesKinData(g.trial).modifiedLocationsInterp{i}(end,3,:)); % ind at which paw reaches max height
                    var{i} = g.sesKinData(g.trial).modifiedLocationsInterp{i}(end,1,maxInd);
                end
            end
            
        case 'stepOverLength'  % length of step over obstacle
            var = getTrialStepLengthsAndDurations(g.sesKinData(g.trial), 0);
        
        case 'preStepOverLength'  % length of step before step over obstacle
            var = getTrialStepLengthsAndDurations(g.sesKinData(g.trial), -1);
            
        case 'prePreStepOverLength'  % length of step before step before step over obstacle
            var = getTrialStepLengthsAndDurations(g.sesKinData(g.trial), -2);
            
        case 'controlStepLength'  % length of final control step (final step before first modified step)
            var = num2cell(nan(1,4));
            if g.sesKinData(g.trial).isTrialAnalyzed
                var = num2cell(g.sesKinData(g.trial).controlSwingLengths(end,:));
            end
        
        case 'isVentralContact'  % whether each paw contacts obs ventrally for >= touchThresh frames
            var = num2cell(nan(1,4));
            if g.sesKinData(g.trial).isTrialAnalyzed
                foreVentInd = find(strcmp(g.sesData.touchClassNames, 'fore_ventral'));
                hindVentInd = find(strcmp(g.sesData.touchClassNames, 'hind_ventral_low'));
                trialTouchesPerPaw = g.sesData.touchesPerPaw(g.sesKinData(g.trial).trialInds,:);
            
                var{1} = sum(trialTouchesPerPaw(:,1)==hindVentInd) >= g.touchThresh;
                var{2} = sum(trialTouchesPerPaw(:,2)==foreVentInd) >= g.touchThresh;
                var{3} = sum(trialTouchesPerPaw(:,3)==foreVentInd) >= g.touchThresh;
                var{4} = sum(trialTouchesPerPaw(:,4)==hindVentInd) >= g.touchThresh;
            end
            
        case 'isDorsalContact'  % whether each paw contacts obs dorsally for >= touchThresh frames
            var = num2cell(nan(1,4));
            if g.sesKinData(g.trial).isTrialAnalyzed
                foreVentInd = find(strcmp(g.sesData.touchClassNames, 'fore_dorsal'));
                hindVentInd = find(strcmp(g.sesData.touchClassNames, 'hind_dorsal'));
                trialTouchesPerPaw = g.sesData.touchesPerPaw(g.sesKinData(g.trial).trialInds,:);
            
                var{1} = sum(trialTouchesPerPaw(:,1)==hindVentInd) >= g.touchThresh;
                var{2} = sum(trialTouchesPerPaw(:,2)==foreVentInd) >= g.touchThresh;
                var{3} = sum(trialTouchesPerPaw(:,3)==foreVentInd) >= g.touchThresh;
                var{4} = sum(trialTouchesPerPaw(:,4)==hindVentInd) >= g.touchThresh;
            end
            
        case 'numTouchFrames'  % number of frames in which paw is in contact with obstacles
            var = num2cell(sum(g.sesData.touchesPerPaw(g.sesKinData(g.trial).trialInds,:)~=0,1));
            
        case 'stepType'  % 1: leading fore // 2: trailing fore // 3: leading hind // 4: traiing hind
            stepTypeSequence = [true false true false;   % leading/trailing
                                true true false false];  % fore/hind
            isLeadingIsLagging = [[g.expData(mouse).sessions(session).trials(trial).paws.isLeading]; ...
                                  [g.expData(mouse).sessions(session).trials(trial).paws.isFore]];  % 2x4 matrix where each column is a paw, and each row is [isLeading, isFore]
            [~, var] = ismember(isLeadingIsLagging', stepTypeSequence', 'rows');
            var = num2cell(var);
            
        case 'distanceToObs'  % linear distance of paw to top of obstacle at zenith of paw trajectory
            var = num2cell(nan(1,4));
            if g.sesKinData(g.trial).isTrialAnalyzed
                for i = find(g.isValidZ(g.trial,:))
                    [~, maxInd] = max(g.sesKinData(g.trial).modifiedLocationsInterp{i}(end,3,:)); % ind at which paw reaches max height
                    xz_paw = g.sesKinData(g.trial).modifiedLocationsInterp{i}(end,[1 3],maxInd);
                    xz_obs = [0 g.sesData.obsHeights(g.trial)/1000];
                    var{i} = norm(xz_paw - xz_obs);  % distance beteen paw and top of obstacle
                end
            end
            
        case 'anyTouchFrames'  % whether the paw contacted the obstacle for any frames
            var = num2cell(any(g.sesData.touchesPerPaw(g.sesKinData(g.trial).trialInds,:)~=0,1));
    end
end



function col = getTableColumn(colName)
    % given column name in sessionInfo, finds values in column name associated with all sessions belonging to a particular mouse
    
    [~, inds] = intersect(g.sessionInfo.session, g.sessions, 'stable');
    col = g.sessionInfo.(colName)(inds);
    if isnumeric(col); col = cellstr(num2str(col)); end
end



function [lengths, durations] = getTrialStepLengthsAndDurations(trialData, stepInd)
    % given kinData for a trial (one row of kinData), returns the lengths 
    % and durations of the steps for all paws for a given stepInd, where 0 
    % is step over, -1 is step before step over, -2 is step before that...
    [lengths, durations] = deal(num2cell(nan(1,4)));
    if trialData.isTrialAnalyzed    
        for i = 1:4
            numModSteps = size(trialData.modifiedLocations{i},1);
            ind = numModSteps+stepInd;
            if ind>0  % if desired step is in modSteps, and not controlSteps
                lengths{i} = trialData.modifiedSwingLengths(ind,i);
                durations{i} = trialData.modifiedSwingDurations(ind,i);
            else
                ind = ind + g.numControlSteps;  % ind of step within controlSteps
                lengths{i} = trialData.controlSwingLengths(ind,i);
                durations{i} = trialData.controlSwingDurations(ind,i);
            end 
        end
    end
end



function avgs = avgSignalPerTrial(sig, sigTimes)
    % averages a signal (e.g. velocity, body angle, tail height) for each
    % trial, ommitting periods after wheel breaks
    
    avgs = cell(1, length(g.sesKinData));
    for i = 1:length(avgs)
        if any(g.sesSpikeData.breaks.times>g.sesData.obsOnTimes(i) & ...
               g.sesSpikeData.breaks.times<g.sesData.obsOffTimes(i))
            endTime = g.sesSpikeData.breaks.times(find(g.sesSpikeData.breaks.times>g.sesData.obsOnTimes(i),1,'first'));
        else
            endTime = g.sesData.obsOffTimes(i);
        end
        avgs{i} = nanmean(sig(sigTimes>g.sesData.obsOnTimes(i) & sigTimes<endTime));
    end
end



function isValidZ = getSessionValidZ(kinData, obsHgts)
    % for each paw in each trial for a session, checks whether paws passes
    % under rather than over the obstacle // returns [trials X paw] logical
    % matrix
    isValidZ = false(length(kinData), 4);
    for i = g.sesKinInds
        for j = 1:4
            atObsInd = find(kinData(i).modifiedLocationsInterp{j}(end,1,:)<=0, 1, 'last'); % final ind at which step is at the obstacle in the x dimension
            if ~isempty(atObsInd)
                isValidZ(i,j) = kinData(i).modifiedLocationsInterp{j}(end,3,atObsInd) > (obsHgts(i)-g.clearanceBuffer); % .001 (m) is a forgiveness factor
            end
        end
    end
end



end
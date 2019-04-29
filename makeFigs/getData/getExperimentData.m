function expData = getExperimentData(sessionInfo, vars, oldData)

% creates nested struct containing data for all mice, sessions, trials, and
% paws // at each level in this heirarchy vars can be computed // add
% sanity checks for all dvs (eg step over hgt cant be less than hgt of asdf
% obstacle...) // dft fcns for getting session, mouse, trial vars, etc...

% TO DO: throw error when var missing, or when var is requested that depends on higher up var that has not been requested...
% make this with recursive function?
% ignore vars if can't be computed (eg if looking for field in sessionInfo that doest exist for given experiment)
% parallelize...

% 'g' stores global variables, allowing things to be passed around!


% settings
metadata = {'touchThresh', 'speedTime', 'preObsLim', 'clearanceBuffer', 'velVsPositionX', 'velContinuousAtContactX'}; % these parameters will be stored as experiment metadata
g.touchThresh = 5;
g.speedTime = .01; % (s) compute velocity over this interval
g.preObsLim = .008; % (m) compute paw height this many meters before it reaches obstacle x postion
g.clearanceBuffer = .001; % (m) trials are excluded in which paw height is less that obsHeight - pawClearnceBuffer at the moment it reaches the x position of the obstacle

g.velVsPositionPrePost = [-.8 .4]; % (m) positions before and after whisker contact to collect velocity data
g.velVsPositionRes = 500; % (tics) how many samples in the x grid

g.velContinuousAtContactPrePost = [-1 1]; % (s) how many seconds before and after wisk contact to compute velocity
g.velContinuousAtContactRes = 500; % (tics) how many samples in the x grid


% initialiations
if ischar(sessionInfo) % if sessionInfo is a string, then it contains the name of a single session, and we will get sesionInfo for this session automatically
    sessionName = sessionInfo;
    sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'sessionInfo.xlsx'), 'Sheet', 'sessions');
    sessionInfo = sessionInfo(strcmp(sessionInfo.session, sessionName),:);
end
g.sessionInfo = sessionInfo;
g.velVsPositionX = linspace(g.velVsPositionPrePost(1), g.velVsPositionPrePost(2), g.velVsPositionRes);
g.velContinuousAtContactX = linspace(g.velContinuousAtContactPrePost(1), g.velContinuousAtContactPrePost(2), g.velContinuousAtContactRes);


mouseVars = {};
sessionVars = {'experiment', 'condition', 'side', 'brainRegion', 'mW', 'conditionNum', 'sessionNum', 'whiskers'};
trialVars = {'obsOnPositions', 'obsOffPositions', 'velContinuousAtContact', 'velVsPosition', 'isLightOn', 'isWheelBreak', 'obsHgt', ...
             'isTrialSuccess', 'trialVel', 'velAtWiskContact', ...
             'trialAngle', 'trialAngleContra', 'angleAtWiskContact', 'angleAtWiskContactContra', ...
             'wiskContactPosition', 'wiskContactTimes', 'isContraFirst', 'isBigStep', 'isModPawContra', ...
             'tailHgt', 'modPawDistanceToObs', 'modPawPredictedDistanceToObs', 'velContinuousAtContact', ...
             'modPawKin', 'modPawKinInterp', 'preModPawKin', 'preModPawKinInterp', 'modPawDeltaLength', 'preModPawDeltaLength', ...
             'sensoryCondition', 'modPawContactInd', 'trialDuration', 'optoOnTimes', 'isOptoOn', 'touchFrames'};
pawVars = {'isContra', 'isFore', 'isLeading', 'isPawSuccess', 'stepOverMaxHgt', 'preObsHgt', 'baselineStepHgt', ...
           'stepOverStartingDistance', 'stepOverEndingDistance', 'stepOverKinInterp', 'controlStepKinInterp', ...
           'isValidZ', 'preObsKin', 'xDistanceAtPeak', 'stepOverLength', 'preStepOverLength', 'prePreStepOverLength'};

% compute only requested vars
if isequal(vars, 'all'); vars = cat(2, sessionVars, trialVars, pawVars); end
mouseVars = mouseVars(ismember(mouseVars, vars));
sessionVars = sessionVars(ismember(sessionVars, vars));
trialVars = trialVars(ismember(trialVars, vars));
pawVars = pawVars(ismember(pawVars, vars));



% get mouse data
g.mice = unique(g.sessionInfo.mouse(logical(g.sessionInfo.include)));
g.expData = struct('mouse', g.mice);
for mouseVar = 1:length(mouseVars)
    temp = getVar(mouseVars{mouseVar}, g);
    [g.expData(1:length(g.mice)).(mouseVars{mouseVar})] = temp{:};
end

% loop over mice
for mouse = 1:length(g.mice)
    g.mouse = mouse;
    
    % get session data
    g.sessions = g.sessionInfo.session(strcmp(g.sessionInfo.mouse, g.mice{mouse}) & logical(g.sessionInfo.include));
    g.expData(mouse).sessions = struct('session', g.sessions);
    for sessionVar = 1:length(sessionVars)
        try
        temp = getVar(sessionVars{sessionVar}, g);
        [g.expData(mouse).sessions(1:length(g.sessions)).(sessionVars{sessionVar})] = temp{:};
        catch; end
    end
    
    % loop over sessions
    for session = 1:length(g.sessions)
        g.session = session;
        
        % check if session already exists in oldData
        if exist('oldData', 'var')
            mouseBin = strcmp({oldData.mouse}, g.mice{mouse});
            if any(mouseBin)
                sesBin = strcmp({oldData(mouseBin).sessions.session}, g.sessions{session});
            else
                sesBin = false;
            end
        else
            sesBin = false;
        end
        
        % if session exists in old data, copy it over
        if any(sesBin)
            fprintf('%s: copying trials from previous data...\n', g.sessions{session})
            g.expData(mouse).sessions(session).trials = oldData(mouseBin).sessions(sesBin).trials;
        
        % otherwise, compute de novo
        else
            % load data from session
            if exist(fullfile(getenv('OBSDATADIR'), 'sessions', g.sessions{session}, 'kinData.mat'), 'file')
                g.sesKinData = load(fullfile(getenv('OBSDATADIR'), 'sessions', g.sessions{session}, 'kinData.mat'), 'kinData');
            else
                try
                    getKinematicData5(g.sessions{session});
                catch
                    fprintf('%s: problem with getKinematicData!', g.sessions{session});
                end
                g.sesKinData = load(fullfile(getenv('OBSDATADIR'), 'sessions', g.sessions{session}, 'kinData.mat'), 'kinData');
            end
            g.sesKinData = g.sesKinData.kinData;
            g.sesKinInds = find([g.sesKinData.isTrialAnalyzed]);
            g.sesData = load(fullfile(getenv('OBSDATADIR'), 'sessions', g.sessions{session}, 'runAnalyzed.mat'));
            warning('off', 'MATLAB:load:variableNotFound');
            g.sesSpikeData = load(fullfile(getenv('OBSDATADIR'), 'sessions', g.sessions{session}, 'run.mat'), 'breaks', 'stimulus');
            g.wheelVel = getVelocity(g.sesData.wheelPositions, g.speedTime, 1000);
            g.isValidZ = getSessionValidZ(g.sesKinData, g.sesData.obsHeightsVid/1000, g);

            % get size of kin data entries
            try
                g.locationsSmps = size(g.sesKinData(find([g.sesKinData.isTrialAnalyzed],1,'first')).modifiedLocations{1}, 3);
                g.locationsInterpSmps = size(g.sesKinData(find([g.sesKinData.isTrialAnalyzed],1,'first')).modifiedLocationsInterp{1}, 3);
            end

            % get trial data
            g.expData(mouse).sessions(session).trials = struct('trial', num2cell(1:length(g.sesKinData)));
            for trialVar = 1:length(trialVars)
                try
                temp = getVar(trialVars{trialVar}, g);
                [g.expData(mouse).sessions(session).trials(1:length(g.sesKinData)).(trialVars{trialVar})] = temp{:};
                catch; end
            end

            % loop over trials
            isStepOverAnalyzed = false(length(g.sesKinData),4);
            for trial = 1:length(g.sesKinData)
                g.trial = trial;

                % get paw data
                g.expData(mouse).sessions(session).trials(trial).paws = struct('paw', {1,2,3,4}, 'pawName', {'LH', 'LF', 'RF', 'RH'});
                for pawVar = 1:length(pawVars)
                    try
                    temp = getVar(pawVars{pawVar}, g);
                    [g.expData(mouse).sessions(session).trials(trial).paws(1:4).(pawVars{pawVar})] = temp{:};
                    catch; end
                end

                if ismember('stepOverKinInterp', vars)
                    isStepOverAnalyzed(trial,:) = cellfun(@(x) ~all(isnan(x(:))), {g.expData(mouse).sessions(session).trials(trial).paws.stepOverKinInterp});
                end
            end

            stepOverAnalyzedRates = num2cell(mean(isStepOverAnalyzed,1));
            fprintf('%s: analysis success rates -> paw1: %.2f, paw2: %.2f, paw3: %.2f, paw4: %.2f \n', ...
                g.sessions{session}, stepOverAnalyzedRates{:})
        end
    end
end
expData = g.expData;
% disp('all done getting experiment data! woo hoo!!!')


% save experiment metadata
dataTemp = expData; clear expData; % nest expData within itself
expData.data = dataTemp; clear dataTemp;
for m = 1:length(metadata); expData.(metadata{m}) = g.(metadata{m}); end % add one field per metadatum









% ---------
% FUNCTIONS
% ---------

function var = getVar(dvName, g) % sessionInfo, expData, mice, mouse, sessions, session, trial, sesKinData, sesData, wheelVel, isValidZ, breakTimes, touchThresh, sesKinInds, locationsInterpSmps, preObsLim
    
    switch dvName
        
        % session variables
        % -----------------
        case 'experiment'
            var = getTableColumn('experiment', g);
        
        case 'condition'
            var = getTableColumn('condition', g);
            
        case 'side'
            var = getTableColumn('side', g);
            
        case 'brainRegion'
            var = getTableColumn('brainRegion', g);
            
        case 'mW'
            var = getTableColumn('mW', g);
            
        case 'conditionNum'
            sessionInfoSub = g.sessionInfo(strcmp(g.sessionInfo.mouse, g.expData(g.mouse).mouse),:);
            var = num2cell(nan(1,height(sessionInfoSub)));
            conditions = unique(sessionInfoSub.condition);
            for i = conditions'
%                 if strcmp(i{1}, 'postContra'); keyboard; end
                bins = strcmp(sessionInfoSub.condition, i{1});
                var(bins) = num2cell(1:sum(bins));
            end
            var = var(logical(sessionInfoSub.include));
            
            
        case 'sessionNum'
            allMouseSessions = g.sessionInfo.session(strcmp(g.sessionInfo.mouse, g.mice{g.mouse})); % all sessions, including those where .include==false
            [~, inds] = intersect(allMouseSessions, g.sessions, 'stable');
            var = num2cell(inds);
            
        case 'whiskers'
            var = getTableColumn('whiskers', g);
        
            
            
           
            
        % trial variables
        % ---------------
        
        case 'obsOnPositions'
            % position of the obstacle relative to mouse nose at the moment it turns on
            var = num2cell(interp1(g.sesData.obsTimes, g.sesData.obsPositionsFixed, g.sesData.obsOnTimes, 'linear'));
            
        case 'obsOffPositions'
            % position of the obstacle relative to mouse nose at the momentit turns off
            var = num2cell(interp1(g.sesData.obsTimes, g.sesData.obsPositionsFixed, g.sesData.obsOffTimes, 'linear'));
        
        case 'velContinuousAtContact'
            % continuous velocity vector surrouding moment of whisker contact
            times = g.velContinuousAtContactX;
            var = repmat({nan(1, length(times))},1,length(g.sesKinData));
            
            for i = g.sesKinInds
                bins = g.sesData.wheelTimes>g.sesData.obsOnTimes(i)-range(times) & ...
                    g.sesData.wheelTimes<g.sesData.obsOffTimes(i)+range(times); % add range(times) as a buffer in case the desired times fall outside the range of the obstacle being on
                var{i} = interp1(g.sesData.wheelTimes(bins), g.wheelVel(bins), times+g.sesData.wiskContactTimes(i));
            end
            
        case 'velVsPosition'
            % continuous velocity as a function of position of obstacle relative to nose
            var = repmat({nan(1, length(g.velVsPositionX))},1,length(g.sesKinData));
            
            for i = g.sesKinInds
                obsAtNoseTime = g.sesData.obsTimes(find(g.sesData.obsPositionsFixed>=0 & g.sesData.obsTimes>g.sesData.obsOnTimes(i), 1, 'first'));
                obsAtNosePos = g.sesData.wheelPositions(find(g.sesData.wheelTimes>obsAtNoseTime,1,'first'));
                trialBins = (g.sesData.wheelPositions > obsAtNosePos+g.velVsPositionPrePost(1)) & (g.sesData.wheelPositions < obsAtNosePos+g.velVsPositionPrePost(2));
                trialPos = g.sesData.wheelPositions(trialBins) - obsAtNosePos; % normalize s.t. 0 corresponds to the position at which the obstacle is at the mouse's nose
                trialVel = g.wheelVel(trialBins); % wheel vel for trial

                % remove duplicate positional values (would be better to average all values in a particular bin)
                [trialPos, uniqueInds] = unique(trialPos, 'stable');
                trialVel = trialVel(uniqueInds);
                
                % interpolate velocities across positional grid and store results
                var{i} = interp1(trialPos, trialVel, g.velVsPositionX, 'linear');
            end
            
        case 'isLightOn'
            var = num2cell(g.sesData.isLightOn);
            
        case 'isWheelBreak'
            var = cell(1,length(g.sesKinData));
            for i = 1:length(g.sesKinData)
               var{i} = any(g.sesSpikeData.breaks.times>g.sesData.obsOnTimes(i) & ...
                            g.sesSpikeData.breaks.times<g.sesData.obsOffTimes(i));
            end
            
        case 'obsHgt'
            var = num2cell(g.sesData.obsHeightsVid/1000); % convert back to meters
        
        case 'isTrialSuccess'
            var = cell(1,length(g.sesKinData));
            for i = 1:length(g.sesKinData)
                var{i} = sum(any(g.sesData.touchesPerPaw(g.sesKinData(i).trialInds,:),2)) < g.touchThresh;
            end
            
        case 'trialVel'
            var = avgSignalPerTrial(g.wheelVel, g.sesData.wheelTimes, g);
        
        case 'velAtWiskContact'
            var = num2cell(interp1(g.sesData.wheelTimes, g.wheelVel, g.sesData.wiskContactTimes));
            
        case 'trialAngle'
            var = avgSignalPerTrial(g.sesData.bodyAngles, g.sesData.frameTimeStamps, g);
            
        case 'trialAngleContra'
            side = g.expData(g.mouse).sessions(g.session).side;
            if strcmp(side, 'left') || strcmp(side, 'right')
                var = getVar('trialAngle', g);
                if strcmp(g.expData(g.mouse).sessions(g.session).side, 'left'); var = num2cell(cellfun(@(x) -x, var)); end % if side is left, then contra limbs are on the right side
            end
        
        case 'angleAtWiskContact'
            bins = ~isnan(g.sesData.frameTimeStamps);
            var = num2cell(interp1(g.sesData.frameTimeStamps(bins), g.sesData.bodyAngles(bins), g.sesData.wiskContactTimes));
            
        case 'angleAtWiskContactContra'
            side = g.expData(g.mouse).sessions(g.session).side;
            if strcmp(side, 'left') || strcmp(side, 'right')
                var = getVar('angleAtWiskContact', g);
                if strcmp(g.expData(g.mouse).sessions(g.session).side, 'left'); var = num2cell(cellfun(@(x) -x, var)); end % if side is left, then contra limbs are on the right side
            end
            
        case 'wiskContactPosition'
            var = num2cell([g.sesKinData.wiskContactPositions]);
            
        case 'wiskContactTimes'
            var = num2cell([g.sesKinData.wiskContactTimes]);
            
        case 'isContraFirst'
            side = g.expData(g.mouse).sessions(g.session).side;
            if strcmp(side, 'left') || strcmp(side, 'right')
                var = num2cell(false(1,length(g.sesKinData)));
                for i = g.sesKinInds; var{i} = ismember(g.sesKinData(i).pawOverSequence(1), [1, 2]); end % is first paw over on left side of body
                if strcmp(side, 'left'); var = num2cell(cellfun(@not, var)); end % if side is left, then contra limbs are on the right side
            end
            
        case 'isBigStep'
            var = num2cell(false(1,length(g.sesKinData)));
            for i = g.sesKinInds
                var{i} = max(g.sesKinData(i).modifiedStepIdentities(:,g.sesKinData(i).firstModPaw))==1; % should i just check the length of modifiedLocations instead?
            end
            
        case 'isModPawContra'
            side = g.expData(g.mouse).sessions(g.session).side;
            if strcmp(side, 'left') || strcmp(side, 'right')
                var = num2cell(false(1,length(g.sesKinData)));
                if strcmp(side, 'right'); contraPaw = 2; else; contraPaw = 3; end
                var([g.sesKinData.isTrialAnalyzed]) = num2cell([g.sesKinData.firstModPaw]==contraPaw);
            end
            
        case 'tailHgt'
            tailHgts = nan(size(g.sesData.frameTimeStamps));
            for i = g.sesKinInds
                tailHgts(g.sesKinData(i).trialInds) = g.sesKinData(i).locationsTail(:,3,1);
            end
            var = avgSignalPerTrial(tailHgts, g.sesData.frameTimeStamps, g);
            
        case 'modPawDistanceToObs'
            var = num2cell(false(1,length(g.sesKinData)));
            for i = g.sesKinInds
                var{i} = g.sesKinData(i).modifiedLocations{g.sesKinData(i).firstModPaw}(1,1,end);
            end
            
        case 'modPawPredictedDistanceToObs'
            var = num2cell(false(1,length(g.sesKinData)));
            for i = g.sesKinInds
                xStart = g.sesKinData(i).modifiedLocations{g.sesKinData(i).firstModPaw}(1,1,1); % first x val of first mod step for mod paw
                predictedLength = g.sesKinData(i).modPredictedLengths(1,g.sesKinData(i).firstModPaw);
                var{i} = xStart + predictedLength;
            end
            
        case 'modPawKin'
            % kinematics of first modified step for first modified paw
            var = repmat({nan(3,g.locationsSmps)},1,length(g.sesKinData));
            
            for i = g.sesKinInds
                var{i} = squeeze(g.sesKinData(i).modifiedLocations{g.sesKinData(i).firstModPaw}(1,:,:));
                if ~g.isValidZ(i,g.sesKinData(i).firstModPaw) && size(g.sesKinData(i).modifiedLocations{g.sesKinData(i).firstModPaw},1)==1 % if step doesn't pass over obstacle and wasn't supposed to pass over obstacle
                    var{i} = nan(3,g.locationsSmps);
                end
            end
            
        case 'modPawKinInterp'
            % kinematics (interpolated) of first modified step for first modified paw
            var = repmat({nan(3,g.locationsInterpSmps)},1,length(g.sesKinData));
            
            for i = g.sesKinInds
                var{i} = squeeze(g.sesKinData(i).modifiedLocationsInterp{g.sesKinData(i).firstModPaw}(1,:,:));
                if ~g.isValidZ(i,g.sesKinData(i).firstModPaw) && size(g.sesKinData(i).modifiedLocations{g.sesKinData(i).firstModPaw},1)==1 % if step doesn't pass over obstacle and wasn't supposed to pass over obstacle
                    var{i} = nan(3,g.locationsInterpSmps);
                end
            end
            
        case 'preModPawKin'
            % kinematics of steps! preceding first modified step for first modified paw
            % note: this matrix contains MULTIPLE steps, and has 3 dimensions as a results ([step X dimension (xyz) X time])
            
            numControlSteps = size(g.sesKinData(find(g.sesKinInds,1,'first')).controlLocations{1},1);
            var = repmat({nan(numControlSteps,3,g.locationsSmps)}, 1, length(g.sesKinData));
            for i = g.sesKinInds
                var{i} = squeeze(g.sesKinData(i).controlLocations{g.sesKinData(i).firstModPaw});
            end
            
        case 'preModPawKinInterp'
            % kinematics (interpolated) of step preceding first modified step for first modified paw
            var = repmat({nan(3,g.locationsInterpSmps)},1,length(g.sesKinData));
            for i = g.sesKinInds
                var{i} = squeeze(g.sesKinData(i).controlLocationsInterp{g.sesKinData(i).firstModPaw}(end,:,:));
            end
            
        case 'modPawDeltaLength'
            % change in length of first modified step of first modified paw
            var = num2cell(false(1,length(g.sesKinData)));
            for i = g.sesKinInds
                predictedLength = g.sesKinData(i).modPredictedLengths(1,g.sesKinData(i).firstModPaw);
                actualLength = g.sesKinData(i).modifiedSwingLengths(1,g.sesKinData(i).firstModPaw);
                var{i} = actualLength - predictedLength;
            end
            
        case 'preModPawDeltaLength' % change in length of step preceding first modified step of first modified paw (serves as control for modPawDeltaLength)
            var = num2cell(false(1,length(g.sesKinData)));
            for i = g.sesKinInds
                predictedLength = g.sesKinData(i).controlPredictedLengths(end,g.sesKinData(i).firstModPaw);
                actualLength = g.sesKinData(i).controlSwingLengths(end,g.sesKinData(i).firstModPaw);
                var{i} = actualLength - predictedLength;
            end
            
        case 'sensoryCondition' % encodes whether animal has light + whiskers ('LW'), whiskers only ('W'), light only ('L'), or neither ('-')
            if ismember('whiskers', g.sessionInfo.Properties.VariableNames)
                hasWhiskers = ~strcmp(g.sessionInfo.whiskers{strcmp(g.sessionInfo.session, g.sessions{g.session})}, 'none');
            else
                hasWhiskers = true;
            end
            if hasWhiskers; conditions = {'W', 'WL'}; else; conditions = {'-', 'L'}; end
            var = conditions([g.sesData.isLightOn]+1);
            
        case 'isValidZ'
            % checks whether step over obstacle passes under obstacle,
            % indicating a tracking error // other variables that rely on z
            % should check whether isValidZ
            
            var = num2cell(g.isValidZ(g.trial,:));
            
        case 'modPawContactInd' % ind in modPawKin at which whiskers contact obstacle
            var = num2cell(nan(1,length(g.sesKinData)));
            for i = g.sesKinInds
                var{i} = g.sesKinData(i).pawObsPosInd(g.sesKinData(i).firstModPaw); 
            end
            
        case 'trialDuration' % duration of trial (from obstalce on to obstacle off time)
            var = num2cell(g.sesData.obsOffTimes - g.sesData.obsOnTimes);
            
        case 'optoOnTimes' % for optogenetics experiments, encodes whether optogenetic stimulation occured for every trial
            light = zscore(g.sesSpikeData.stimulus.values);
            times = g.sesSpikeData.stimulus.times;
            
            var = num2cell(nan(1,length(g.sesKinData)));
            for i = 1:length(g.sesKinData)
%                 trialLight = light(times>g.sesData.obsOnTimes(i) & times<g.sesData.obsOffTimes(i));
%                 var{i} = (sum(trialLight>1)*g.sesSpikeData.stimulus.interval) > .01; % if there is greater than 10 ms of light
                lightOnInd = find(times>g.sesData.obsOnTimes(i) & times<g.sesData.obsOffTimes(i) & light>1, 1, 'first');
                if ~isempty(lightOnInd); var{i} = times(lightOnInd); end
            end
            
        case 'isOptoOn'
            var = num2cell(~cellfun(@isnan, getVar('optoOnTimes', g)));
            
        case 'touchFrames'
            var = cell(1,length(g.sesKinData));
            for i = 1:length(g.sesKinData)
                var{i} = sum(any(g.sesData.touchesPerPaw(g.sesKinData(i).trialInds,:),2));
            end
            
            
            

            
        % paw variables
        % -------------
        case 'isContra'
            side = g.expData(g.mouse).sessions(g.session).side;
            if strcmp(side, 'left') || strcmp(side, 'right')
                var = [1 1 0 0]; % left side contra by default
                if strcmp(g.expData(g.mouse).sessions(g.session).side, 'left'); var = ~var; end % if side is left, then contra limbs are on the right side
                var = num2cell(var);
            end
            
        case 'isFore'
            var = num2cell(logical([0 1 1 0]));
        
        case 'isLeading'
            var = num2cell(logical([1 1 0 0])); % start assuming left is leading
            seq = g.sesKinData(g.trial).pawOverSequence;
            if find(seq==4)<find(seq==1); var([1,4]) = var([4,1]); end
            if find(seq==3)<find(seq==2); var([2,3]) = var([3,2]); end
            
        case 'isPawSuccess'
            var = num2cell(sum(g.sesData.touchesPerPaw(g.sesKinData(g.trial).trialInds,:),1) < g.touchThresh);
        
        case 'stepOverMaxHgt'
            % maximum height of paw on step over obstacle
            
            var = num2cell(nan(1,4));
            if g.sesKinData(g.trial).isTrialAnalyzed
                var(g.isValidZ(g.trial,:)) = num2cell(cellfun(@(x) max(x(end,3,:)), ...
                    g.sesKinData(g.trial).modifiedLocationsInterp(g.isValidZ(g.trial,:))));
                if any(cellfun(@(x) x<0, var(g.isValidZ(g.trial,:)))); keyboard; end
            end
            
        case 'preObsHgt' % height of paw when the step over is preObsLim in front of osbtacle
            var = num2cell(nan(1,4));
            if g.sesKinData(g.trial).isTrialAnalyzed
                for i = find(g.isValidZ(g.trial,:))
                    ind = find(g.sesKinData(g.trial).modifiedLocationsInterp{i}(end,1,:)>-g.preObsLim, 1, 'first');
                    if ind>1; var{i} = g.sesKinData(g.trial).modifiedLocationsInterp{i}(end,3,ind); end
                end
            end
            
        case 'baselineStepHgt'
            var = num2cell(nan(1,4));
            if g.sesKinData(g.trial).isTrialAnalyzed
                var = cellfun(@(x) nanmean(max(squeeze(x(:,3,:)),[],2)), ... % avg the max z for each noObsStep
                    g.sesKinData(g.trial).noObsLocations, 'UniformOutput', false);
            end
            
%         case 'penultStepLength'
%             var = num2cell(nan(1,4));
%             if g.sesKinData(g.trial).isTrialAnalyzed
%                 for i = 1:4
%                     if max(g.sesKinData(g.trial).modifiedStepIdentities(:,i))==1 % if 1 mod step, penultimate step is final control step
%                         var{i} = g.sesKinData(g.trial).controlSwingLengths(end,i);
%                     else
%                         var{i} = g.sesKinData(g.trial).modifiedSwingLengths(end-1,i); % otherwise it is the second to last mod step
%                     end 
%                 end
%             end
            
        case 'stepOverStartingDistance' % distance of paw to obs at start of step over obs
            var = num2cell(nan(1,4));
            if g.sesKinData(g.trial).isTrialAnalyzed
                var = num2cell(cellfun(@(x) x(end,1,1), g.sesKinData(g.trial).modifiedLocations));
            end
            
        case 'stepOverEndingDistance' % distance of paw to obs at end of step over obs
            var = num2cell(nan(1,4));
            if g.sesKinData(g.trial).isTrialAnalyzed
                var = num2cell(cellfun(@(x) x(end,1,end), g.sesKinData(g.trial).modifiedLocationsInterp));
            end
        
        case 'stepOverKinInterp' % interpolated kinematics of step over obstacle
            var = repmat({nan(3,g.locationsInterpSmps)},1,4);
            if g.sesKinData(g.trial).isTrialAnalyzed
                var(g.isValidZ(g.trial,:)) = cellfun(@(x) squeeze(x(end,:,:)), ...
                    g.sesKinData(g.trial).modifiedLocationsInterp(g.isValidZ(g.trial,:)), 'UniformOutput', false);
            end
            
        case 'controlStepKinInterp' % interpolated kinematics of step over obstacle
            var = repmat({nan(3,g.locationsInterpSmps)},1,4);
            if g.sesKinData(g.trial).isTrialAnalyzed
                var(g.isValidZ(g.trial,:)) = cellfun(@(x) squeeze(x(end,:,:)), ...
                    g.sesKinData(g.trial).controlLocationsInterp(g.isValidZ(g.trial,:)), 'UniformOutput', false);
            end
            
        case 'preObsKin' % interpolated kinematics of step over obs, but only before step reaches obstacle
            var = repmat({nan(3,g.locationsInterpSmps)},1,4);
            if g.sesKinData(g.trial).isTrialAnalyzed
                for i = find(g.isValidZ(g.trial,:))
                    ind = find(g.sesKinData(g.trial).modifiedLocations{i}(end,1,:)>-g.preObsLim, 1, 'first');
                    kin = squeeze(g.sesKinData(g.trial).modifiedLocations{i}(end,:,1:ind));
                    
                    if ind>1
                        for j = 1:3 % x, y, z
                            var{i}(j,:) = interp1(1:size(kin,2), kin(j,:), linspace(1,size(kin,2), g.locationsInterpSmps), 'linear');
                        end 
                    end
                end
            end
            
        case 'xDistanceAtPeak' % distance of paw to obstacle when it is at max height
            var = num2cell(nan(1,4));
            if g.sesKinData(g.trial).isTrialAnalyzed
                for i = find(g.isValidZ(g.trial,:))
                    [~, maxInd] = max(g.sesKinData(g.trial).modifiedLocationsInterp{i}(end,3,:)); % ind at which paw reaches max height
                    var{i} = g.sesKinData(g.trial).modifiedLocationsInterp{i}(end,1,maxInd);
                end
            end
            
        case 'stepOverLength' % length of step over obstacle
            var = getTrialStepLengthsAndDurations(g.sesKinData(g.trial), 0);
        
        case 'preStepOverLength' % length of step over obstacle
            var = getTrialStepLengthsAndDurations(g.sesKinData(g.trial), -1);
            
        case 'prePreStepOverLength' % length of step over obstacle
            var = getTrialStepLengthsAndDurations(g.sesKinData(g.trial), -2);
    end
end




function col = getTableColumn(colName, g)
    % given column name in sessionInfo, finds values in column name associated with all sessions belonging to a particular mouse
    
    [~, inds] = intersect(g.sessionInfo.session, g.sessions, 'stable');
    col = g.sessionInfo.(colName)(inds);
    if isnumeric(col); col = cellstr(num2str(col)); end
end



function [lengths, durations] = getTrialStepLengthsAndDurations(trialData, stepInd)
    % give kinData for a trial (one row of kinData), returns the lengths of
    % the steps for all paws for a given stepInd, where 0 is step over, -1
    % is step before step over, -2 is step before that, etc...
    [lengths, durations] = deal(num2cell(nan(1,4)));
    numControlSteps = size(trialData.controlSwingDurations,1);
    if trialData.isTrialAnalyzed    
        for i = 1:4
            numModSteps = size(trialData.modifiedLocations{i},1);
            ind = numModSteps+stepInd;
            if ind>0 % if desired step is in modSteps, and not controlSteps
                lengths{i} = trialData.modifiedSwingLengths(ind,i);
                durations{i} = trialData.modifiedSwingDurations(ind,i);
            else
                ind = numControlSteps + (numModSteps + stepInd); % ind of step within controlSteps
                lengths{i} = trialData.controlSwingLengths(ind,i);
                durations{i} = trialData.controlSwingDurations(ind,i);
            end 
        end
    end
end



function avgs = avgSignalPerTrial(sig, sigTimes, g)
    % averages a signal (velocity, body angle, tail height) for each trial, but ommitting periods after wheel breaks for each trial
    
    avgs = cell(1, length(g.sesKinData));
    for i = 1:length(avgs)
        
        if any(g.sesSpikeData.breaks.times>g.sesData.obsOnTimes(i) & ...
               g.sesSpikeData.breaks.times<g.sesData.obsOffTimes(i))
            endTime = g.sesSpikeData.breaks.times(find(g.sesSpikeData.breaks.times>g.sesData.obsOnTimes(i),1,'first'));
        else
            endTime = g.sesData.obsOffTimes(i);
        end
        
        inds = sigTimes>g.sesData.obsOnTimes(i) & sigTimes<endTime;
        avgs{i} = nanmean(sig(inds));
    end
end



function isValidZ = getSessionValidZ(kinData, obsHgts, g)
    % for each paw in each trial for a session, checks whether paws passes
    % under rather than over the obstacle // returns [trials X paw] logical
    % matrix
    
    isValidZ = false(length(kinData), 4);
    for i = g.sesKinInds
        for j = 1:4
            atObsInd = find(kinData(i).modifiedLocationsInterp{j}(end,1,:)>=0, 1, 'first'); % ind at which step is at the obstacle in the x dimension
            if ~isempty(atObsInd)
                isValidZ(i,j) = kinData(i).modifiedLocationsInterp{j}(end,3,atObsInd) > (obsHgts(i)-g.clearanceBuffer); % .001 (m) is a forgiveness factor
            end
        end
    end
end




end
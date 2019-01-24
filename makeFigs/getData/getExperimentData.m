function expData = getExperimentData(sessionInfo, vars)

% creates nested struct containing data for all mice, sessions, trials, and
% paws // at each level in this heirarchy vars can be computed


% settings
touchThresh = 5;
speedTime = .05; % compute velocity over this interval


% initialiations
mouseVars = {};
sessionVars = {'condition', 'side', 'brainRegion'};
trialVars = {'isLightOn', 'isWheelBreak', 'obsHgt', 'isTrialSuccess', 'velAtWiskContact'}; % , 'velAvg', 'angleAvg', 'angleAtContact', 'trialVel', 'isContraFirst', 'isBigStep', 'wiskContactPos', 'tailHeight', 'obsPosAtContact', 'modPawStepNum'
pawVars = {'stepOverMaxHeight'}; % 'penultLength', 'isPawSuccess', 'stepOverKin', 'preObsHeight', 'baselineHgt', 'firstModKin'

% compute only requested vars
mouseVars = mouseVars(ismember(mouseVars, vars));
sessionVars = sessionVars(ismember(sessionVars, vars));
trialVars = trialVars(ismember(trialVars, vars));
pawVars = pawVars(ismember(pawVars, vars));


% get mouse data
mice = unique(sessionInfo.mouse(logical(sessionInfo.include)));
expData = struct('mouse', mice);
for mouseVar = 1:length(mouseVars)
    temp = getVar(mouseVars{mouseVar});
    [expData(1:length(mice)).(mouseVars{mouseVar})] = temp{:};
end

% loop over mice
for mouse = 1:length(mice)
    
    % get session data
    sessions = sessionInfo.session(strcmp(sessionInfo.mouse, mice{mouse}));
    expData(mouse).sessions = struct('session', sessions);
    for sessionVar = 1:length(sessionVars)
        temp = getVar(sessionVars{sessionVar});
        [expData(mouse).sessions(1:length(sessions)).(sessionVars{sessionVar})] = temp{:};
    end
    
    % loop over sessions
    for session = 1:length(sessions)
        
        % load data from session
        fprintf('%s: computing\n', sessions{session})
        sesKinData = load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{session}, 'kinData.mat'), 'kinData'); sesKinData = sesKinData.kinData;
        sesData = load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{session}, 'runAnalyzed.mat'));
        wheelVel = getVelocity(sesData.wheelPositions, speedTime, 1000);
        
        % get trial data
        expData(mouse).sessions(session).trials = struct('trial', num2cell(1:length(sesKinData)));
        for trialVar = 1:length(trialVars)
            temp = getVar(trialVars{trialVar});
            [expData(mouse).sessions(session).trials(1:length(sesKinData)).(trialVars{trialVar})] = temp{:};
        end
        
        % loop over trials
        for trial = 1:length(sesKinData)
            
            % get paw data
            expData(mouse).sessions(session).trials(trial).paws = struct('paw', {1,2,3,4}, 'pawName', {'LH', 'LF', 'RF', 'RH'});
            for pawVar = 1:length(pawVars)
                temp = getVar(pawVars{pawVar});
                [expData(mouse).sessions(session).trials(trial).paws(1:4).(pawVars{pawVar})] = temp{:};
            end
        end
    end 
end
disp('all done getting experiment data! woo hoo!!!')









% ---------
% FUNCTIONS
% ---------

function var = getVar(dvName)
    
    switch dvName
        
        % session variables
        % -----------------
        case 'condition'
            var = getTableColumn('condition');
            
        case 'side'
            var = getTableColumn('side');
            
        case 'brainRegion'
            var = getTableColumn('brainRegion');
            
        % trial variables
        % ---------------
        case 'isLightOn'
            var = num2cell(sesData.isLightOn);
            
        case 'isWheelBreak'
            load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{session}, 'run.mat'), 'breaks');
            var = cell(1,length(sesKinData));
            for i = 1:length(sesKinData)
               var{i} = any(breaks.times>sesData.obsOnTimes(i) & breaks.times<sesData.obsOffTimes(i));
            end
            
        case 'obsHgt'
            var = num2cell(sesData.obsHeightsVid);
        
        case 'isTrialSuccess'
            var = cell(1,length(sesKinData));
            for i = 1:length(sesKinData)
                var{i} = sum(any(sesData.touchesPerPaw(sesKinData(i).trialInds,:),2)) < touchThresh;
            end
            
        case 'velAtWiskContact'
            var = num2cell(interp1(sesData.wheelTimes, wheelVel, sesData.wiskContactTimes));
            
        % paw variables
        % -------------
        case 'stepOverMaxHeight'
            if ~isempty(sesKinData(trial).modifiedLocations)
                var = num2cell(cellfun(@(x) max(x(end,3,:)), sesKinData(trial).modifiedLocations));
            else
               var = num2cell(nan(1,4));
            end
    end
end




function col = getTableColumn(colName)
    % given column name in sessionInfo, finds values in column name associated with all sessions belonging to a particular mouse
    [~, inds] = intersect(sessionInfo.session, sessions, 'stable');
    col = sessionInfo.(colName)(inds);
end
end
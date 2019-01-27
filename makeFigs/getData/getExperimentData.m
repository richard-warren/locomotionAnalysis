function expData = getExperimentData(sessionInfo, vars)

% creates nested struct containing data for all mice, sessions, trials, and
% paws // at each level in this heirarchy vars can be computed // add
% sanity checks for all dvs (eg step over hgt cant be less than hgt of asdf
% obstacle...) // dft fcns for getting session, mouse, trial vars, etc...

% TO DO: throw error when var missing, or when var is requested that depends on higher up var that has not been requested... // make this with recursive function? // ignore vars if can't be computed (eg if
% looking for field in sessionInfo that doest exist for given experiment)

% settings
touchThresh = 5;
speedTime = .05; % compute velocity over this interval


% initialiations
mouseVars = {};
sessionVars = {'condition', 'side', 'brainRegion'};
trialVars = {'isLightOn', 'isWheelBreak', 'obsHgt', 'isTrialSuccess', 'trialVel', 'velAtWiskContact', ...
             'trialAngle', 'trialAngleContra', 'angleAtWiskContact', 'angleAtWiskContactContra', ...
             'obsPosAtContact', 'wiskContactPositions', 'isContraFirst', 'isBigStep'};
pawVars = {'isContra', 'pawType', 'isPawSuccess', 'stepOverMaxHgt', 'baselineStepHgt', 'penultStepLength'};

% 'preObsHeight', 'stepOverKin', 'firstModKin'

% compute only requested vars
if isequal(vars, 'all'); vars = cat(2, sessionVars, trialVars, pawVars); end
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
        sesKinInds = find([sesKinData.isTrialAnalyzed]);
        sesData = load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{session}, 'runAnalyzed.mat'));
        sesBreaks = load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{session}, 'run.mat'), 'breaks'); breakTimes = sesBreaks.breaks.times;
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
            var = cell(1,length(sesKinData));
            for i = 1:length(sesKinData)
               var{i} = any(breakTimes>sesData.obsOnTimes(i) & breakTimes<sesData.obsOffTimes(i));
            end
            
        case 'obsHgt'
            var = num2cell(sesData.obsHeightsVid);
        
        case 'isTrialSuccess'
            var = cell(1,length(sesKinData));
            for i = 1:length(sesKinData)
                var{i} = sum(any(sesData.touchesPerPaw(sesKinData(i).trialInds,:),2)) < touchThresh;
            end
            
        case 'trialVel'
            var = avgSignalPerTrial(wheelVel, sesData.wheelTimes);
        
        case 'velAtWiskContact'
            var = num2cell(interp1(sesData.wheelTimes, wheelVel, sesData.wiskContactTimes));
            
        case 'trialAngle'
            var = avgSignalPerTrial(sesData.bodyAngles, sesData.frameTimeStamps);
            
        case 'trialAngleContra'
            var = getVar('trialAngle');
            if strcmp(expData(mouse).sessions(session).side, 'left'); var = num2cell(cellfun(@(x) -x, var)); end % if side is left, then contra limbs are on the right side
        
        case 'angleAtWiskContact'
            bins = ~isnan(sesData.frameTimeStamps);
            var = num2cell(interp1(sesData.frameTimeStamps(bins), sesData.bodyAngles(bins), sesData.wiskContactTimes));
            
        case 'angleAtWiskContactContra'
            var = getVar('angleAtWiskContact');
            if strcmp(expData(mouse).sessions(session).side, 'left'); var = num2cell(cellfun(@(x) -x, var)); end % if side is left, then contra limbs are on the right side
            
        case 'obsPosAtContact'
            var = num2cell(interp1(sesData.obsTimes, sesData.obsPositionsFixed, sesData.wiskContactTimes));
            
        case 'wiskContactPositions'
            var = num2cell([sesKinData.wiskContactPositions]);
            
        case 'isContraFirst'
            var = num2cell(false(1,length(sesKinData)));
            for i = sesKinInds; var{i} = ismember(sesKinData(i).pawOverSequence(1), [1, 2]); end % is first paw over on left side of body
            if strcmp(expData(mouse).sessions(session).side, 'left'); var = num2cell(cellfun(@not, var)); end % if side is left, then contra limbs are on the right side
            
        case 'isBigStep'
            var = num2cell(false(1,length(sesKinData)));
            for i = sesKinInds
                var{i} = max(sesKinData(i).modifiedStepIdentities(:,sesKinData(i).firstModPaw))==1; %
            end
            
        
        
        % paw variables
        % -------------
        case 'isContra'
            var = [1 1 0 0]; % left side contra by default
            if strcmp(expData(mouse).sessions(session).side, 'left'); var = ~var; end % if side is left, then contra limbs are on the right side
            var = num2cell(var);
        
        case 'pawType'
            var = cell(1,4);
            seq = sesKinData(trial).pawOverSequence;
            if find(seq==1)<find(seq==4); var([1,4]) = {'leadHind', 'lagHind'}; else; var([1,4]) = {'lagHind', 'leadHind'}; end
            if find(seq==2)<find(seq==3); var([2,3]) = {'leadFore', 'lagFore'}; else; var([2,3]) = {'lagFore', 'leadFore'}; end
            
        case 'isPawSuccess'
            var = num2cell(sum(sesData.touchesPerPaw(sesKinData(trial).trialInds,:),1) < touchThresh);
        
        case 'stepOverMaxHgt'
            var = num2cell(nan(1,4));
            if sesKinData(trial).isTrialAnalyzed
                var = num2cell(cellfun(@(x) max(x(end,3,:)), sesKinData(trial).modifiedLocations));
            end
            
        case 'baselineStepHgt'
            var = num2cell(nan(1,4));
            if sesKinData(trial).isTrialAnalyzed
                var = cellfun(@(x) nanmean(max(squeeze(x(:,3,:)),[],2)), ... % avg the max z for each noObsStep
                    sesKinData(trial).noObsLocations, 'UniformOutput', false);
            end
            
        case 'penultStepLength'
            var = num2cell(nan(1,4));
            if sesKinData(trial).isTrialAnalyzed
                for i = 1:4
                    if max(sesKinData(trial).modifiedStepIdentities(:,i))==1 % if 1 mod step, penultimate step is final control step
                        var{i} = sesKinData(trial).controlSwingLengths(end,i);
                    else
                        var{i} = sesKinData(trial).modifiedSwingLengths(end-1,i); % otherwise it is the second to last mod step
                    end 
                end
            end
    end
end




function col = getTableColumn(colName)
    % given column name in sessionInfo, finds values in column name associated with all sessions belonging to a particular mouse
    
    [~, inds] = intersect(sessionInfo.session, sessions, 'stable');
    col = sessionInfo.(colName)(inds);
end



function avgs = avgSignalPerTrial(sig, sigTimes)
    % averages a signal (velocity, body angle, tail height) for each trial, but ommitting periods after wheel breaks for each trial
    
    avgs = cell(1, length(sesKinData));
    for i = 1:length(avgs)
        
        if any(breakTimes>sesData.obsOnTimes(i) & breakTimes<sesData.obsOffTimes(i))
            endTime = breakTimes(find(breakTimes>sesData.obsOnTimes(i),1,'first'));
        else
            endTime = sesData.obsOffTimes(i);
        end
        
        inds = sigTimes>sesData.obsOnTimes(i) & sigTimes<endTime;
        avgs{i} = nanmean(sig(inds));
    end
end




end
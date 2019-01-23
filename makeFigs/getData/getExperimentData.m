function data = getExperimentData(vars, sessionInfo, sessionData)

% creates nested struct containing data for all mice, sessions, trials, and
% paws // at each level in this heirarchy vars can be computed

% temp
% vars = {'condition', 'isTrialSuccess', 'baselineHgt', 'isPawSuccess'};
% getDv = @(x) 1;


% settings
touchThresh = 5;


% initialiations
mouseVars = {};
sessionVars = {'condition', 'side', 'brainRegion'};
trialVars = {'isLightOn', 'isTrialSuccess', 'vel', 'bodyAngle', 'isContraFirst', ... % !!! vel
    'isBigStep', 'wiskContactPos', 'obsHgt', 'tailHeight'};
pawVars = {'maxHeight', 'penultLength', 'isPawSuccess', 'stepOverKinenatics', 'preObsHeight', 'baselineHgt'};

% compute only requested vars
mouseVars = mouseVars(ismember(mouseVars, vars));
sessionVars = sessionVars(ismember(sessionVars, vars));
trialVars = trialVars(ismember(trialVars, vars));
pawVars = pawVars(ismember(pawVars, vars));


% get mouse data
mice = unique(sessionInfo.mouse(logical(sessionInfo.include)));
data = struct('mouse', mice);
for mouseVar = 1:length(mouseVars)
    temp = getVar(mouseVars{mouseVar});
    [data(1:length(mice)).(mouseVars{mouseVar})] = temp{:};
end

% loop over mice
for mouse = 1:length(mice)
    
    % get session data
    sessions = sessionInfo.session(strcmp(sessionInfo.mouse, mice{mouse}));
    data(mouse).sessions = struct('session', sessions);
    % !!! load session data here
    for sessionVar = 1:length(sessionVars)
        temp = getVar(sessionVars{sessionVar});
        [data(mouse).sessions(1:length(sessions)).(sessionVars{sessionVar})] = temp{:};
    end
    
    % loop over sessions
    for session = 1:length(sessions)
        
        % get trial data
        data(mouse).sessions(session).trials = struct('trial', num2cell(1:length(sessionData)));
        for trialVar = 1:length(trialVars)
            temp = getVar(trialVars{trialVar});
            [data(mouse).sessions(session).trials(1:length(sessionData)).(trialVars{trialVar})] = temp{:};
        end
        
        % loop over trials
        for trial = 1:length(sessionData)
            
            % get paw data
            data(mouse).sessions(session).trials(trial).paws = struct('paw', {1,2,3,4}, 'pawName', {'LH', 'LF', 'RF', 'RH'});
            for pawVar = 1:length(pawVars)
                temp = getVar(pawVars{pawVar});
                data(mouse).sessions(session).trials(trial).paws(1:length(paws)).(pawVars{pawVar}) = temp{:};
            end
            
        end
    end 
end




function var = getVar(dvName)
    
    switch dvName
        
        % session variables
        % -----------------
        case 'condition'
            [~, inds] = intersect(sessionInfo.session, sessions, 'stable');
            var = sessionInfo.condition(inds);
            
        case 'side'
            [~, inds] = intersect(sessionInfo.session, sessions, 'stable');
            var = sessionInfo.side(inds);
            
        case 'brainRegion'
            [~, inds] = intersect(sessionInfo.session, sessions, 'stable');
            var = sessionInfo.brainRegion(inds);
            
        % trial variables
        % ---------------
        case 'isLightOn'
            var = num2cell([sessionData.isLightOn]);
        
        case 'isTrialSuccess'
            var = cell(1,length(sessionData));
            for i = 1:length(sessionData); var{i} = sum(any(sessionData(i).trialTouchesPerPaw,2)) < touchThresh; end
            
        % paw variables
        % -------------
            
    end
end
end
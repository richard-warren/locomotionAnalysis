function sessionDvs = getSessionDvs(dvs, speedAvoidanceData, kinData)

% given speedAvoidanceData and kinData structures (which contain trial to
% trial information for multiple sessions), computes dependent measures
% averaged across each session // only computes dependent measures listed
% in dvs cell array // can also include desired metadata fieds (eg
% 'conditionNum'), and they will be added to sessionDvs

% TO DO: add minTrial criterion // make this work when only speedAvoidanceData is given // add ability to avg using weights to create matched distributions


% settings
touchThresh = 5; % if paw contacts obs for more than touchThresh frames, trial is considered touching

% initializations
kinSessions = unique({kinData.session});
speedAvoidanceSessions = unique({speedAvoidanceData.session});
if isequal(kinSessions, speedAvoidanceSessions); sessions = kinSessions; else; disp('WARNING! different sessions detected in kinData and speedAvoidanceData!'); end


% get avg for each dv for each session
for i = 1:length(sessions)
    
    sessionBinsKin = strcmp({kinData.session}, sessions{i});
    sessionBinsSpdAv = strcmp({speedAvoidanceData.session}, sessions{i});
    
    % store session metadata
    sessionDvs(i).session = sessions{i};
    sessionDvs(i).mouse = speedAvoidanceData(find(sessionBinsSpdAv,1,'first')).mouse;
    sessionDvs(i).condition = speedAvoidanceData(find(sessionBinsSpdAv,1,'first')).condition;
    sessionDvs(i).conditionNum = speedAvoidanceData(find(sessionBinsSpdAv,1,'first')).conditionNum;
    sessionDvs(i).sessionNum= speedAvoidanceData(find(sessionBinsSpdAv,1,'first')).sessionNum;
    
    % store dependent variables
    for j = 1:length(dvs)
        sessionDvs(i).(dvs{j}) = getDv(dvs{j});
    end
end


% function that returns a single dependent variable per call
function dv = getDv(dvName)
    
    switch dvName
        
        case 'success'
            isSuccess = cellfun(@(x) sum(any(x,2)), {speedAvoidanceData(sessionBinsSpdAv).trialTouchesPerPaw}) < touchThresh;
            dv = nanmean(isSuccess);
            
        case 'speed'
            dv = nanmean([speedAvoidanceData(sessionBinsSpdAv).avgVel]);
            
        case 'body_angle'
            dv = nanmean([speedAvoidanceData(sessionBinsSpdAv).avgAngle]);
            side = unique({speedAvoidanceData(sessionBinsSpdAv).side});
            if strcmp(side, 'left'); dv = -dv; end
        
        case 'conditionNum'
            dv = unique([speedAvoidanceData(sessionBinsSpdAv).conditionNum]);
            
        case 'ipsiErrRate'
            side = unique({speedAvoidanceData(sessionBinsSpdAv).side});
            leftErrorRate = nanmean(cellfun(@(x) sum(any(x(:,[1 2]),2))>=touchThresh, {speedAvoidanceData(sessionBinsSpdAv).trialTouchesPerPaw}));
            rightErrorRate = nanmean(cellfun(@(x) sum(any(x(:,[3 4]),2))>=touchThresh, {speedAvoidanceData(sessionBinsSpdAv).trialTouchesPerPaw}));
            if strcmp(side, 'left'); dv = leftErrorRate; else; dv = rightErrorRate; end
        
        case 'contraErrRate'
            side = unique({speedAvoidanceData(sessionBinsSpdAv).side});
            leftErrorRate = nanmean(cellfun(@(x) sum(any(x(:,[1 2]),2))>=touchThresh, {speedAvoidanceData(sessionBinsSpdAv).trialTouchesPerPaw}));
            rightErrorRate = nanmean(cellfun(@(x) sum(any(x(:,[3 4]),2))>=touchThresh, {speedAvoidanceData(sessionBinsSpdAv).trialTouchesPerPaw}));
            if strcmp(side, 'left'); dv = rightErrorRate; else; dv = leftErrorRate; end
    end
end



end
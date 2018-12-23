function sessionDvs = getSessionDvs(dvs, speedAvoidanceData, kinData)

% given speedAvoidanceData and kinData structures (which contain trial to
% trial information for multiple sessions), computes dependent measures
% averaged across each session // only computes dependent measures listed
% in dvs cell array // can also include desired metadata fieds (eg
% 'conditionNum'), and they will be added to sessionDvs // kinData is only
% necessary if some of the dvs dependent upon data in this structure

% TO DO: add minTrial criterion // make this work when only speedAvoidanceData is given // add ability to avg using weights to create matched distributions


% settings
touchThresh = 5; % if paw contacts obs for more than touchThresh frames, trial is considered touching

% initializations
if exist('kinData', 'var')
    kinSessions = unique({kinData.session});
    speedAvoidanceSessions = unique({speedAvoidanceData.session});
    if isequal(kinSessions, speedAvoidanceSessions); sessions = kinSessions; else; disp('WARNING! different sessions detected in kinData and speedAvoidanceData!'); end
else
    sessions = unique({speedAvoidanceData.session});
end

% get avg for each dv for each session
for i = 1:length(sessions)
    
    if exist('kinData', 'var'); sessionBinsKin = strcmp({kinData.session}, sessions{i}); end
    sessionBinsSpdAv = strcmp({speedAvoidanceData.session}, sessions{i});
    if isfield(speedAvoidanceData, 'side'); side = unique({speedAvoidanceData(sessionBinsSpdAv).side}); end % many dvs depend on knowing which side of brain/whiskers were manipulated // compute this once hear to avoid rdundant computation downstream
    
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
            
        case 'bodyAngleContra'
            dv = nanmean([speedAvoidanceData(sessionBinsSpdAv).avgAngle]);
            if strcmp(side, 'left'); dv = -dv; end
        
        case 'conditionNum'
            dv = unique([speedAvoidanceData(sessionBinsSpdAv).conditionNum]);
            
        case 'ipsiErrRate'
            leftErrorRate = nanmean(cellfun(@(x) sum(any(x(:,[1 2]),2))>=touchThresh, {speedAvoidanceData(sessionBinsSpdAv).trialTouchesPerPaw}));
            rightErrorRate = nanmean(cellfun(@(x) sum(any(x(:,[3 4]),2))>=touchThresh, {speedAvoidanceData(sessionBinsSpdAv).trialTouchesPerPaw}));
            if strcmp(side, 'left'); dv = leftErrorRate; else; dv = rightErrorRate; end
        
        case 'contraErrRate'
            leftErrorRate = nanmean(cellfun(@(x) sum(any(x(:,[1 2]),2))>=touchThresh, {speedAvoidanceData(sessionBinsSpdAv).trialTouchesPerPaw}));
            rightErrorRate = nanmean(cellfun(@(x) sum(any(x(:,[3 4]),2))>=touchThresh, {speedAvoidanceData(sessionBinsSpdAv).trialTouchesPerPaw}));
            if strcmp(side, 'left'); dv = rightErrorRate; else; dv = leftErrorRate; end
            
        case 'contraFirstRate'
            dv = mean([kinData(sessionBinsKin).firstPawOver]==2);
            if strcmp(side, 'left'); dv = 1-dv; end
    end
end



end
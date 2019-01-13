function sessionDvs = getSessionDvs(dvs, speedAvoidanceData, kinData)

% given speedAvoidanceData and kinData structures (which contain trial to
% trial information for multiple sessions), computes dependent measures
% averaged across each session // only computes dependent measures listed
% in dvs cell array // can also include desired metadata fieds (eg
% 'conditionNum'), and they will be added to sessionDvs // kinData is only
% necessary if some of the dvs dependent upon data in this structure

% many dvs compute the dv for all trials, then broken down by ipsi and
% contra trials (e.g. leading limb paw height overall, and broken down by
% ipsi/contra)

% TO DO: add ability to avg using weights to create matched distributions // add minTrial criterion


% settings
touchThresh = 5; % if paw contacts obs for more than touchThresh frames, trial is considered touching
preObsLim = .008; % determine paw height when it is preObsLim meters in from of obstacle

% initializations
if exist('kinData', 'var')
    kinSessions = unique({kinData.session});
    speedAvoidanceSessions = unique({speedAvoidanceData.session});
    if isequal(kinSessions, speedAvoidanceSessions); sessions = kinSessions; else; disp('WARNING! different sessions detected in kinData and speedAvoidanceData!'); end
else
    sessions = unique({speedAvoidanceData.session});
end

% get info about paw heights only if a dv needs it
if any(ismember(dvs, {'pawHgt', 'hgtShaping'}))
    pawHgts = nan(1,length(kinData));
    for i = 1:length(kinData)
        firstPawOver = kinData(i).firstPawOver;
        locations = squeeze(kinData(i).modifiedLocations{firstPawOver}(end,:,:));
        hgtInd = find(locations(1,:)>-preObsLim,1,'first');
        if hgtInd>1; pawHgts(i) = locations(3, hgtInd)*1000; end
    end
    obsHgts = [kinData.obsHeightsVid];
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

% create ipsi and contra fields for all dvs with overall, ipsi, and contra versions
sessionDvs = struct2table(sessionDvs);
for i = 1:length(dvs)
    if size(sessionDvs.(dvs{i}),2)==3  % if a dv has an overall, ipsi, and contra version
        try
        sessionDvs.([dvs{i} 'Ipsi']) = sessionDvs.(dvs{i})(:,2);
        sessionDvs.([dvs{i} 'Contra']) = sessionDvs.(dvs{i})(:,3);
        sessionDvs.(dvs{i}) = sessionDvs.(dvs{i})(:,1);
        catch; keyboard; end
    end
end
sessionDvs = table2struct(sessionDvs);



% --------
% function that returns a single dependent variable per call
% --------

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
            
        case 'forePawErrRate'
            errorRate = nanmean(cellfun(@(x) sum(any(x(:,[2 3]),2))>=touchThresh, {speedAvoidanceData(sessionBinsSpdAv).trialTouchesPerPaw}));
            leftErrorRate = nanmean(cellfun(@(x) sum(any(x(:,2),2))>=touchThresh, {speedAvoidanceData(sessionBinsSpdAv).trialTouchesPerPaw}));
            rightErrorRate = nanmean(cellfun(@(x) sum(any(x(:,3),2))>=touchThresh, {speedAvoidanceData(sessionBinsSpdAv).trialTouchesPerPaw}));
            dv = [errorRate, leftErrorRate, rightErrorRate]; % [overal, ipsi, contra]
            if strcmp(side, 'right'); dv = dv([1 3 2]); end
            
        case 'contraFirstRate'
            dv = mean([kinData(sessionBinsKin).firstPawOver]==2);
            if strcmp(side, 'left'); dv = 1-dv; end
            
        case 'bigStepProb'
            modStepNum = cellfun(@(stepNums, paw) stepNums(paw), ...
                {kinData(sessionBinsKin).modStepNum}, {kinData(sessionBinsKin).firstModPaw});
            probOverall = mean(modStepNum==1);
            probLeft = mean(modStepNum([kinData(sessionBinsKin).firstModPaw]==2)==1);
            probRight = mean(modStepNum([kinData(sessionBinsKin).firstModPaw]==3)==1);
            dv = [probOverall, probLeft, probRight];
            if strcmp(side, 'right'); dv = dv([1 3 2]); end
            
        case 'pawHgt'
            hgtOverall = nanmean(pawHgts(sessionBinsKin));
            hgtLeft = nanmean(pawHgts(sessionBinsKin & [kinData.firstModPaw]==2));
            hgtRight = nanmean(pawHgts(sessionBinsKin & [kinData.firstModPaw]==3));
            dv = [hgtOverall, hgtLeft, hgtRight];
            if strcmp(side, 'right'); dv = dv([1 3 2]); end
            
        case 'hgtShaping'
            nonNans = ~isnan(obsHgts) & ~isnan(pawHgts);
            bins = {sessionBinsKin & nonNans, ...
                    sessionBinsKin & nonNans & [kinData.firstModPaw]==2, ...
                    sessionBinsKin & nonNans & [kinData.firstModPaw]==3}; % bins for overall, left, and right
            
            polyFits = cell(1,3);
            for k = 1:3
                polyFits{k} = polyfit(obsHgts(bins{k}), pawHgts(bins{k}), 1);
            end
            dv = [polyFits{1}(1), polyFits{2}(1), polyFits{3}(1)];
            if strcmp(side, 'right'); dv = dv([1 3 2]); end
            
    end
end



end
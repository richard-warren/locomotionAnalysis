% function getObsTrajectories(sessions)

% temp
sessions = {'180122_001', '180122_002', '180122_003'};%, ...
%             '180123_001', '180123_002', '180123_003', ...
%             '180124_001', '180124_002', '180124_003'}; ...
%             '180125_001', '180125_002', '180125_003'};

% settings
obsPos = -.0073;

% initializations
sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx']);
data = struct();
dataInd = 1;

for i = 1:length(sessions)
    
    % report progress
    fprintf('%s: collecting data\n', sessions{i});
    
    % load session data
    load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\runAnalyzed.mat'],...
            'obsPositions', 'obsTimes', 'obsPixPositions', 'frameTimeStamps', ...
            'obsOnTimes', 'obsOffTimes', 'nosePos', 'targetFs', 'wheelPositions', 'wheelTimes');
    obsPositions = fixObsPositions(obsPositions, obsTimes, obsPixPositions, frameTimeStamps, obsOnTimes, obsOffTimes, nosePos(1));
    load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\tracking\locationsBotCorrected.mat'], 'locations')
    load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\tracking\stanceBins.mat'], 'stanceBins')
    
    % get velocities for all trials in session
    sessionVels = getTrialSpeedsAtObsPos(obsPos, wheelPositions, wheelTimes, obsPositions, obsTimes, obsOnTimes, speedTime, targetFs);
    
    for j = 1:length(obsOnTimes)
        
        % get swing identities
        
        
        
        
        % store results
        sessionInfoBin = find(strcmp(sessionInfo.session, sessions{i}));
        data(dataInd).mouse = sessionInfo.mouse{sessionInfoBin};
        data(dataInd).session = sessions{i};
        data(dataInd).vel = sessionVels(j);
        dataInd = dataInd + 1;
        
    end
end

fprintf('--- done collecting data ---\n', sessions{i});










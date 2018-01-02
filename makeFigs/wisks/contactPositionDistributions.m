


% temp
mice = {'run6', 'run7', 'run8'};


% user settings
sessiosnToInclude = 3; % for the plot comparing light vs no light speed, only take the most recent lightVsNoLightSessions for each mouse
wiskTouchThresh = -.75; % if wiskTouchSignal surpasses wiskTouchThresh, then wiskTouchPixels are drawn on wisk frame to show points of contact
obsNosePos = .333;

% initializations
sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx']);

sessionBins = ismember(sessionInfo.mouse, mice) &...
              strcmp(sessionInfo.experiment, 'obsBr') &...
              sessionInfo.include;
sessions = sessionInfo(sessionBins, :);

data = struct(); % stores trial data for all sessions
condColors = flipud(winter(2)); % light on, light off;


%% iterate through sessions

for i = 1:size(sessions,1)
    
    % load session data
    load([getenv('OBSDATADIR') 'sessions\' sessions.session{i} '\runAnalyzed.mat'],...
            'obsPixPositions', 'frameTimeStamps',...
            'obsPositions', 'obsTimes',...
            'obsOnTimes', 'obsOffTimes',...
            'obsLightOnTimes', 'obsLightOffTimes',...
            'wiskTouchSignal', 'frameTimeStampsWisk');
    obsPositions = fixObsPositions(obsPositions, obsTimes, obsOnTimes); % correct for drift in obstacle position readings
    
    
    % convert wisk contacts to z scores
    realInds = ~isnan(wiskTouchSignal);
    normedReal = zscore(wiskTouchSignal(realInds));
    wiskTouchSignal = nan(size(wiskTouchSignal));
    wiskTouchSignal(realInds) = normedReal;
    
    
    % get pix position of first contact for all trial
    contactPositions = nan(length(obsOnTimes), 1);
    isLightOn = false(length(obsOnTimes), 1);
    
    for j = 1:length(obsOnTimes)
        
        % get position of first contact
        touchTime = frameTimeStampsWisk(find(frameTimeStampsWisk>obsOnTimes(j) & wiskTouchSignal>wiskTouchThresh, 1, 'first'));
        if ~isempty(touchTime)
%             [minDiff, posInd] = min(abs(obsTimes-touchTime));
%             if minDiff < .002
%                 contactPositions(j) = obsPositions(posInd);
%             end
              
            frameInd = find(abs(frameTimeStamps-touchTime)<.002);
            if ~isempty(frameInd)
                contactPositions(j) = obsPixPositions(frameInd);
            end
        end
        
        % find whether light was on
        isLightOn(j) = min(abs(obsOnTimes(j) - obsLightOnTimes)) < 1; % did the light turn on near whether the obstacle turned on
        
    end
    
    data(i).mouse = sessions.mouse{i};
    data(i).lightOnContactPositions = contactPositions(isLightOn);
    data(i).lightOffContactPositions = contactPositions(~isLightOn);   
end



%% plot results
close all; figure('color', [1 1 1]);

% get all contact positions
lightOnContacts = {data.lightOnContactPositions}; lightOnContacts = cat(1, lightOnContacts{:});
lightOffContacts = {data.lightOffContactPositions}; lightOffContacts = cat(1, lightOffContacts{:});

% convert to real world units
load('mToPix', 'mapping'); % load linear mapping to meters to Pixels
lightOnContactsMm = lightOnContacts / mapping(1) * 1000;
lightOffContactsMm = lightOffContacts / mapping(1) * 1000;
meanPos = nanmean([lightOnContactsMm; lightOffContactsMm]);

histogram(lightOnContactsMm-meanPos, 25, 'facecolor', condColors(1,:), 'normalization', 'probability'); hold on
histogram(lightOffContactsMm-meanPos, 25,'facecolor', condColors(2,:), 'normalization', 'probability')

% pimp fig
set(gca, 'box', 'off')










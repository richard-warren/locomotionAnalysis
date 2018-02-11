


% temp
sessions = {'180122_000', '180122_001', '180122_002', '180122_003', ...
            '180123_000', '180123_001', '180123_002', '180123_003', ...
            '180124_000', '180124_001', '180124_002', '180124_003', ...
            '180125_000', '180125_001', '180125_002', '180125_003'};


% settings
wiskTouchThresh = -.75; % if wiskTouchSignal surpasses wiskTouchThresh, then wiskTouchPixels are drawn on wisk frame to show points of contact
obsNosePos = .336; % !!! these scripts will fail if you try to incorporate sessions in which the camera settings are different of the position of the mouse/headplate have changed

% initializations
sessionInfo = readtable([getenv('OBSDATADIR') 'sessions\sessionInfo.xlsx']);
data(length(sessions)) = struct(); % stores trial data for all sessions



%% iterate through sessions

for i = 1:length(sessions)
    
    disp(sessions{i})
    
    % load session data
    load([getenv('OBSDATADIR') 'sessions\' sessions{i} '\runAnalyzed.mat'],...
            'obsPixPositions', 'frameTimeStamps',...
            'obsPositions', 'obsTimes',...
            'obsOnTimes', 'obsOffTimes',...
            'obsLightOnTimes', 'obsLightOffTimes',...
            'wiskTouchSignal', 'frameTimeStampsWisk', 'mToPixMapping');
    vidWisk = VideoReader([getenv('OBSDATADIR') 'sessions\' sessions{i} '\runWisk.mp4']);
    vidTop = VideoReader([getenv('OBSDATADIR') 'sessions\' sessions{i} '\runTop.mp4']);
    obsPositions = fixObsPositions(obsPositions, obsTimes, obsOnTimes); % correct for drift in obstacle position readings
    
    
    % convert wisk contacts to z scores
    realInds = ~isnan(wiskTouchSignal);
    normedReal = zscore(wiskTouchSignal(realInds));
    wiskTouchSignal = nan(size(wiskTouchSignal));
    wiskTouchSignal(realInds) = normedReal;
    
    
    % get pix position of first contact for all trial
    contactPositions = nan(length(obsOnTimes), 1);
    data(i).contactFramesWisk = uint8(nan(vidWisk.Height, vidWisk.Width, length(obsOnTimes)));
    data(i).contactFramesTop = uint8(nan(vidTop.Height, vidTop.Width, length(obsOnTimes)));
    
    for j = 1:length(obsOnTimes)
        
        % get position of first contact
        contactIndWisk = find(frameTimeStampsWisk>obsOnTimes(j) & wiskTouchSignal>wiskTouchThresh, 1, 'first');
        contactTime = frameTimeStampsWisk(contactIndWisk);
        
        if ~isempty(contactTime)
            
            contactIndTop = find(abs(frameTimeStamps-contactTime)<.002);
            
            if ~isempty(contactIndTop)
                
                contactPositions(j) = obsPixPositions(contactIndTop);
                
                contactFrameWisk = rgb2gray(read(vidWisk, contactIndWisk));
                data(i).contactFramesWisk(:,:,j) = contactFrameWisk;
                
                contactFrameTop = rgb2gray(read(vidTop, contactIndTop));
                data(i).contactFramesTop(:,:,j) = contactFrameTop;
                
            end
        end
    end
    
    % convert to meters
    mToPixMapping = median(mToPixMapping,1);
    contactPositions = (contactPositions-mToPixMapping(2)) / mToPixMapping(1);
    
    % subtract nose position
    contactPositions = contactPositions - obsNosePos;
    
    % store data
    sessionInfoBin = find(strcmp(sessionInfo.session, sessions{i}));
    data(i).mouse = sessionInfo.mouse{sessionInfoBin};
    data(i).contactPositions = contactPositions;
    data(i).mToPixMapping = mToPixMapping;
    
end



%% plot results
close all; figure('color', [1 1 1]);

% get all contact positions
allContactPositions = {data.contactPositions};
allContactPositions = cat(1, allContactPositions{:});

% get average contact frame
allContactFrames = {data.contactFramesWisk};
allContactFrames = cat(3, allContactFrames{:});

% convert to real world units
histogram(allContactPositions, 25, 'normalization', 'probability'); hold on

% pimp fig
set(gca, 'box', 'off')

%% show all contact frames












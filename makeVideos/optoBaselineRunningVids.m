% for each session, renders one video per light condition showing running
% around the time of stimulation


% settings
sessions = {'190807_000', '190807_001', '190807_002'};
powers = [0.12, 0.35, 1.0];
trialsPerVid = 5;
targetFps = 50;
timePrePost = [-.5 2]; % time before and after opto onset to show
baseDir = fullfile(getenv('OBSDATADIR'), 'editedVid', 'opto', 'noObstacles');


% initializations
sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', 'optoNotes');
sessionInfo = sessionInfo(ismember(sessionInfo.session, sessions),:); % remove empty rows, not included sessions
mice = unique(sessionInfo.mouse);
brainRegions = unique(sessionInfo.brainRegion);

% collect data all data in big ol table
for i = 1:length(sessions)
    
    % load session data
    fprintf('%s: rendering running opto videos...\n', sessions{i})
    load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'runAnalyzed.mat'), ...
        'wheelPositions', 'wheelTimes', 'obsOnTimes', 'obsOffTimes', 'targetFs', 'frameTimeStamps')
    load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'run.mat'), 'stimulus')
    vidTop = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'runTop.mp4'));
    vidBot = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'runBot.mp4'));
    
    % get light conditions for each trial
    trialPowers = nan(1, length(obsOnTimes));
    for j = 1:length(obsOnTimes)
        trialStim = stimulus.values(stimulus.times>obsOnTimes(j) & stimulus.times<obsOffTimes(j));
        trialPower = max(trialStim)/5;  % peak signal is trialPower fraction of 5V max
        [minDif, minInd] = min(abs(powers-trialPower));  % find closest power in powers, defined above
        if minDif<.01; trialPowers(j) = powers(minInd); end  % trialPower is nan if close value is not in powers, defined above
    end
    sessionPowers = unique(trialPowers);
    sessionPowers = sessionPowers(~isnan(sessionPowers));
    
    % render videos, omg
    for j = 1:length(sessionPowers)
        
        sessionBin = strcmp(sessionInfo.session, sessions{i});
        fileName = fullfile(baseDir, sprintf('%s, %s, %s, %.2fpower.mp4', ...
            sessions{i}, sessionInfo.mouse{sessionBin}, sessionInfo.brainRegion{sessionBin}, sessionPowers(j)));
        vidWriter = VideoWriter(fileName, 'MPEG-4');
        set(vidWriter, 'FrameRate', targetFps);
        open(vidWriter);
        
        trialsToShow = find(trialPowers==sessionPowers(j));
        trialsToShow = sort(trialsToShow(randsample(length(trialsToShow), trialsPerVid)));
        
        for k = trialsToShow 
            trialInds = find(frameTimeStamps>(obsOnTimes(k)+timePrePost(1)) & ...
                             frameTimeStamps<(obsOnTimes(k)+timePrePost(2)));
            for m = trialInds'
                frame = rgb2gray(cat(1, read(vidTop, m), read(vidBot, m)));
                text = sprintf('trial %i', k);
                frame = insertText(frame, [size(frame,2) size(frame,1)], text,...
                                   'BoxColor', 'black', 'AnchorPoint', 'RightBottom', 'TextColor', 'white');
                writeVideo(vidWriter, frame);
            end
        end
        
    end
end








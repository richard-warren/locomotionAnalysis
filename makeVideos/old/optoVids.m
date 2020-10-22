% for each session, renders one video per light condition showing running
% around the time of stimulation

% todo: get fastest trials only?


% settings
sessions = {'190822_000', '190822_001', '190822_002'};  % mtc, 2mm dorsal
% sessions = {'190823_000', '190823_001', '190823_002'};  % sen, 2mm dorsal
% sessions = {'190826_000', '190826_001', '190826_002'};  % olf, 2mm dorsal

trialsPerVid = 15;
targetFps = 50;
timePreOpto = -.25; % time before opto to show (trial continues until obs offset)
baseDir = fullfile(getenv('OBSDATADIR'), 'editedVid', 'opto');


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
    sessionBin = strcmp(sessionInfo.session, sessions{i});
    powers = cellfun(@str2num, strsplit(sessionInfo.power___{sessionBin}, ', '), 'UniformOutput', false);
    powers = [0, cat(2, powers{:})];
    if contains(sessionInfo.experiment{sessionBin}, 'noObs'); folder='noObstacles'; else; folder='obstacles'; end
    
    data = getExperimentData(sessions{i}, 'all');
    sesPowers = powers(knnsearch(powers', [data.data.sessions.trials.optoPower]'));
    optoOnTimes = [data.data.sessions.trials.optoOnTimes];
    
    % render videos, omg
    for j = 1:length(powers)
        
        fileName = fullfile(baseDir, folder, sprintf('%s, %s, %s, %.2fpower.mp4', ...
            sessions{i}, sessionInfo.mouse{sessionBin}, sessionInfo.brainRegion{sessionBin}, powers(j)));
        vidWriter = VideoWriter(fileName, 'MPEG-4');
        set(vidWriter, 'FrameRate', targetFps);
        open(vidWriter);
        
        trialsToShow = find(sesPowers==powers(j));
        trialsToShow = sort(trialsToShow(randsample(length(trialsToShow), min(trialsPerVid, length(trialsToShow)))));
        
        if powers(j)==0; startTimes = obsOnTimes; else;  startTimes = optoOnTimes; end
        
        for k = trialsToShow 
            trialInds = find(frameTimeStamps>(startTimes(k)+timePreOpto) & ...
                             frameTimeStamps<(obsOffTimes(k)));
            for m = trialInds'
                frame = rgb2gray(cat(1, read(vidTop, m), read(vidBot, m)));
                text = sprintf('trial %i', k);
                frame = insertText(frame, [size(frame,2) size(frame,1)], text,...
                    'BoxColor', 'black', 'AnchorPoint', 'RightBottom', 'TextColor', 'white');
                if frameTimeStamps(m)>startTimes(k) && powers(j)>0
                    frame = insertText(frame, [0 size(frame,1)], 'opto on',...
                        'BoxColor', 'white', 'AnchorPoint', 'LeftBottom', 'TextColor', 'black', 'BoxOpacity', 1);
                end
                writeVideo(vidWriter, frame);
            end
        end
        close(vidWriter)
    end
end
disp('all done!')




%% randomize for nate

% settings
folder = fullfile(getenv('OBSDATADIR'), 'editedVid', 'opto', 'obstacles'); % folder containing obstacle vids

% make new folder for these data
dateTime = datestr(datetime(now, 'ConvertFrom', 'datenum'), 'yymmdd-HH.MM.SS');
blindTestFolder = fullfile(folder, 'blindTests', dateTime);
mkdir(blindTestFolder);


vids = dir(fullfile(folder, '*.mp4'));

sessionNames = cellfun(@(x) x(1:10), {vids.name}, 'UniformOutput', false);
sessions = unique(sessionNames);
mouseNames = cellfun(@(x) x(13:16), {vids.name}, 'UniformOutput', false);
mice = unique(mouseNames);

vidIdStruct = struct();  % struct containing randomly assigned session and vid IDs
ind = 1;

for m = 1:length(mice)
    
    mouseSessions = unique(sessionNames(strcmp(mouseNames, mice{m})));
    sessionIds = randsample(length(mouseSessions), length(mouseSessions));  % randomly assigned IDs for each sessions
    
    for i = 1:length(sessionIds)
        
        vidInds = find(contains({vids.name}, mouseSessions{i}));
        vidIds = randsample(length(vidInds), length(vidInds));
        
        for j = 1:length(vidIds)
            vidOriginal = vids(vidInds(j)).name;
            vidNew = sprintf('%s_ses%i_vid%i.mp4', mice{m}, sessionIds(i), vidIds(j));
            copyfile(fullfile(folder, vidOriginal), ...
                     fullfile(blindTestFolder, vidNew));
            
            vidIdStruct(ind).vidNew = vidNew;
            vidIdStruct(ind).vidOriginal = vidOriginal;
            ind = ind + 1;
        end
    end
end

% save vid IDs to spreadsheet
writetable(struct2table(vidIdStruct), ...
    fullfile(blindTestFolder, 'key.csv'));
disp('all done!')








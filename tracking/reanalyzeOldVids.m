%% find out which sessions have 'originalDimensions' files

files = dir(fullfile(getenv('OBSDATADIR'), 'sessions'));
sessions = {files([files.isdir]).name};
origSessions = {};

fileExists = false(1,length(sessions));
fprintf('\n\n--------------------looking for originalDimensions--------------------\n')
for i = 1:length(sessions)
    dirSub = dir(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, '*.mp4'));
    bins = contains({dirSub.name}, 'originalDimensions');
    if any(bins)
        fprintf('%s: ', sessions{i})
        fprintf('%s ', dirSub(bins).name)
        fprintf('\n')
        origSessions{end+1} = sessions{i};
    end
end
disp('all done!')

%% find ephysSessions

ephysInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'ephysInfo.xlsx'));
ephysSessions = ephysInfo.session(ephysInfo.include==1);
clear ephysInfo

% %% find sessions with ephys folders
% 
% unusableSessions = {'181109_000', ...  % ephys folder only
%                     '190923_003', ...  % ephys folder only
%                     '190523_000', ...
%                     '190523_001', ...
%                     '190523_002', ...
%                     '200118_000'};
% 
% files = dir(fullfile(getenv('OBSDATADIR'), 'sessions'));
% sessions = {files([files.isdir]).name};
% ephysSessions = {};
% 
% hasEphysFolder = false(1,length(sessions));
% fprintf('\n\n--------------------looking for ephys folders--------------------\n')
% for i = 1:length(sessions)
%     dirSub = dir(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}));
%     dirSub = dirSub([dirSub.isdir]);
%     bins = contains({dirSub.name}, 'ephys_');
%     if any(bins)
%         fprintf('%s: ', sessions{i})
%         fprintf('%s ', dirSub(bins).name)
%         fprintf('\n')
%         ephysSessions{end+1} = sessions{i};
%     end
% end
% disp('all done!')


% %% get rid of cropped views and concat top and bot FOR EPHYS SESSIONS ONLY
% % note: this does not delete the un-concatenated versions, which could optionally be done later to save disk space
% 
% for i = 1:length(ephysSessions)
%     
%     folder = fullfile(getenv('OBSDATADIR'), 'sessions', ephysSessions{i});
%     dirSub = dir(fullfile(folder, '*.mp4'));
%     
%     % rename originalDimensions
%     origInds = find(contains({dirSub.name}, '_originalDimensions'));
%     if ~isempty(origInds)
%         fprintf('%s: renaming files ', ephysSessions{i})
%         for j = 1:length(origInds)
%             fprintf('%s ', dirSub(origInds(j)).name)
%             movefile(fullfile(folder, dirSub(origInds(j)).name), ...
%                      fullfile(folder, erase(dirSub(origInds(j)).name, '_originalDimensions')));
%         end
%         fprintf('\n')
%         
%         % concatenate views if originally recorded un-concatenated
%         if exist(fullfile(folder, 'runTop.mp4')); concatTopBotVids(ephysSessions{i}); end
%         fprintf('\n')
%     end
% end

%% reanalyze everything for ephys sessions

problemSessions = {};

for i = 1:length(ephysSessions)
    fprintf('\n___________%i/%i___________\n', i, length(ephysSessions))
    
    % concat top and bot if necessary
    folder = fullfile(getenv('OBSDATADIR'), 'sessions', ephysSessions{i});
    if ~exist(fullfile(folder, 'run.mp4')); concatTopBotVids(ephysSessions{i}); end
    
    try
        analyzeSession(ephysSessions{i}, ...
            'overwriteVars', 'all', ...
            'verbose', true, ...
            'superVerbose', false, ...
            'rerunRunNetwork', true, ...
            'rerunWiskNetwork', true, ...
            'rerunPawContactNetwork', true, ...
            'rerunWiskContactNetwork', true);
    catch
        fprintf('%s: problem with analysis!\n', ephysSessions{i})
        problemSessions{end+1} = ephysSessions{i};
    end
end
disp('all done!')


%% reanalyze single field in ephysSessions

% settings
skipSessions = {};
vars = {'lickTimes'};

sessions = ephysSessions(~ismember(ephysSessions, skipSessions));
for i = 1:length(sessions)
    analyzeSession(sessions{i}, 'overwriteVars', vars, 'verbose', true);
    fprintf('\n')
end
disp('all done!')

%% show tracking with continuous signal

session = ephysSessions{1};

locationsWisk = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw_wisk.csv'));
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), 'frameTimeStampsWisk')

showTracking(session, 'sig', locationsWisk.tongue_1, 'sigTimes', frameTimeStampsWisk)


%% reanalyze single sessions
analyzeSession('191008_003', ...
            'overwriteVars', 'all', ...
            'verbose', true, ...
            'superVerbose', false, ...
            'rerunRunNetwork', true, ...
            'rerunWiskNetwork', true, ...
            'rerunPawContactNetwork', true, ...
            'rerunWiskContactNetwork', true);

%% recover broken session

session = '191009_003'; % '200118_001', '191009_003'

% figure out if any frames lost at the beginning of session
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), 'ledInds', 'ledIndsWisk', 'rewardTimes')
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'run.mat'), 'exposure')
firstLedTtl = find(exposure.times > rewardTimes(1), 1, 'first');

if firstLedTtl==ledInds(1) && firstLedTtl==ledIndsWisk(1)
    disp('no frames lost at beginning')
else
    fprintf('%i frames unaccounted for in run camera\n', firstLedTtl-ledInds(1))
    fprintf('%i frames unaccounted for in wisk camera\n', firstLedTtl-ledIndsWisk(1))
end

% check length of videos
vid = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'run.mp4'));
vidWisk = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runWisk.mp4'));
fprintf('%i run frames, %i whisker frames\n', vid.NumberOfFrames, vidWisk.NumberOfFrames)

% check if there are skipped frames
runCamMeta = dlmread(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'run.csv')); % columns: bonsai timestamps, point grey counter, point grey timestamps (uninterpretted)
wiskCamMeta = dlmread(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'wisk.csv')); % columns: bonsai timestamps, point grey counter, point grey timestamps (uninterpretted)
deltaFramesRun = diff(runCamMeta(:,2));
deltaFramesWisk = diff(wiskCamMeta(:,1));

if any(deltaFramesRun>1) || any(deltaFramesWisk>1)
    disp('frames were skipped!')
else
    disp('no skipped frames')
end

% save data assuming no frames lost at beginning
data = load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'));

data.frameTimeStamps = nan(length(runCamMeta),1);
data.frameTimeStampsWisk = nan(length(wiskCamMeta),1);
data.frameTimeStamps(1:length(exposure.times)) = exposure.times;
data.frameTimeStampsWisk(1:length(exposure.times)) = exposure.times;

save(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), '-struct', 'data')
disp('data saved')


%% to clear up disk space we could:


% get rid of runTop, runBot when run exists

% get rid of originalDimensions for old sessions OR ...

% get rid of cropped for old sessions and reanalyze: could potentially
% analyze on cropped vids with old network, then shift the coordinates to
% accomodate the cropping, and throw away cropped vid, performing the rest
% of the analysis on the uncropped video // alternatively, see if old DLC
% can handle uncropped vids, and reanalyze like that...

%% show sample whisker frame for each session

for i = 1:length(ephysSessions)
    vid = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', ephysSessions{i}, 'runWisk.mp4'));
    figure('name', ephysSessions{i}); imshow(read(vid,10000));
end

%% play video while zooming in on face



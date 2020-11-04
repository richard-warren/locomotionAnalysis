%% prepare tracking data for sample session for andrea

% settings
session = '181215_003';
outdir = 'C:\Users\richa\Desktop\adrea_files';

% load tracking data
[pawXYZ, pawXYZ_pixels] = getPawXYZ(session);  % in pixels
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'rewardTimes', 'obsOnTimes', 'obsOffTimes', 'frameTimeStamps')

% pack locations into a matrix
[paws, paws_pixels] = deal(nan(length(frameTimeStamps), 4, 3));
names = {'paw1LH', 'paw2LF', 'paw3RF', 'paw4RH'};
for i = 1:4
    paws(:,i,:) = pawXYZ.(names{i});
    paws_pixels(:,i,:) = pawXYZ_pixels.(names{i});
end

% save files
obstacleTimes = [obsOnTimes, obsOffTimes];
t = frameTimeStamps;
save(fullfile(outdir, ['trackingData_' session '.mat']), ...
    'paws', 'paws_pixels', 'obstacleTimes', 't', 'rewardTimes')

% copy raw tracking data
copyfile(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw.csv'), ...
    fullfile(outdir, 'rawTracking.csv'))

%% prepare tracking data for sample session for andrea

% settings
session = '200901_000';
outdir = 'C:\Users\richa\Desktop\adrea_files';

% load tracking data
sesFolder = fullfile(getenv('OBSDATADIR'), 'sessions', session);
% [pawXYZ, pawXYZ_pixels] = getPawXYZ(session);  % in pixels
load(fullfile(sesFolder, 'runAnalyzed.mat'), ...
    'rewardTimes', 'obsOnTimes', 'obsOffTimes', 'frameTimeStamps', ...
    'wheelPositions', 'wheelTimes')

%% pack locations into a matrix
[paws, paws_pixels] = deal(nan(length(frameTimeStamps), 4, 3));
names = {'paw1LH', 'paw2LF', 'paw3RF', 'paw4RH'};
for i = 1:4
    paws(:,i,:) = pawXYZ.(names{i});
    paws_pixels(:,i,:) = pawXYZ_pixels.(names{i});
end

% save files
obstacleTimes = [obsOnTimes, obsOffTimes];
t = frameTimeStamps;
save(fullfile(outdir, ['trackingData_' session '.mat']), ...
    'paws', 'paws_pixels', 'obstacleTimes', 't', 'rewardTimes')

% copy files
copyFiles = {'run.mp4', 'runWisk.mp4'};
for i = 1:length(copyFiles)
    copyfile(fullfile(sesFolder, copyFiles{i}), fullfile(outdir, copyFiles{i}));
end
%%





















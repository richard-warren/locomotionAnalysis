%% prepare tracking data for sample session for andrea

% settings
sessions = getEphysSessions();
sessions = sessions(end-6:end);


for s = 1:length(sessions)
    session = sessions{s};
    outdir = fullfile('C:\Users\rick\Desktop\andrea_files', session);
    mkdir(outdir)
    fprintf('prepping %s...\n', session)
    
    % load tracking data
    sesFolder = fullfile(getenv('OBSDATADIR'), 'sessions', session);
    [pawXYZ, pawXYZ_pixels] = getPawXYZ(session);  % in pixels
    load(fullfile(sesFolder, 'runAnalyzed.mat'), ...
        'rewardTimes', 'obsOnTimes', 'obsOffTimes', 'frameTimeStamps', ...
        'wheelPositions', 'wheelTimes', 'whiskerAngle', 'frameTimeStampsWisk', ...
        'lickTimes', 'bodyAngles', 'wiskContactTimes')
    locationsWisk = readtable(fullfile(sesFolder, 'trackedFeaturesRaw_wisk.csv'));
    t = frameTimeStamps;
    obstacleTimes = [obsOnTimes, obsOffTimes];

    % pack locations into a matrix
    [paws, paws_pixels] = deal(nan(length(frameTimeStamps), 4, 3));
    names = {'paw1LH', 'paw2LF', 'paw3RF', 'paw4RH'};
    for i = 1:4
        paws(:,i,:) = pawXYZ.(names{i});
        paws_pixels(:,i,:) = pawXYZ_pixels.(names{i});
    end

    % wheel velocity
    vel = getVelocity(wheelPositions, .05, 1/diff(wheelTimes(1:2)));
    vel = interp1(wheelTimes, vel, t);

    % whiskerAngle
    whiskerAngle = interp1(frameTimeStampsWisk, whiskerAngle, t);

    % jaw
    [jawX, jawZ, jawConfidence] = deal(locationsWisk.jaw, locationsWisk.jaw_1, locationsWisk.jaw_2);
    jawX(jawX<prctile(jawX,1)) = nan; jawX(jawX>prctile(jawX,99)) = nan;
    jawZ(jawZ<prctile(jawZ,1)) = nan; jawZ(jawZ>prctile(jawZ,99)) = nan;
    jaw = [jawX jawZ];
    bins = all(~isnan(jaw),2) & jawConfidence>.8;
    jaw(~bins, :) = nan;
    jaw = fillmissing(jaw, 'linear');
    jaw = interp1(frameTimeStampsWisk, jaw, t);

    save(fullfile(outdir, ['trackingData_' session '.mat']), ...
        'paws', 'paws_pixels', 'obstacleTimes', 't', 'rewardTimes', ...
        'vel', 'jaw', 'whiskerAngle', 'bodyAngles', 'wiskContactTimes', ...
        'lickTimes')

    % copy files
    copyFiles = {'run.mp4', 'runWisk.mp4'};
    for i = 1:length(copyFiles)
        copyfile(fullfile(sesFolder, copyFiles{i}), fullfile(outdir, copyFiles{i}));
    end
end
function plotEphysLocomotorScatters3(session, opts)

% settings
s.folder = fullfile(getenv('OBSDATADIR'), 'figures', 'ephys', 'LocomotorScatters', session);  % folder in which the plots will be saved
% s.firingRatePencitiles = [30, 30]; % only include modified steps with the lowest 30% FR and highest 30% FR
s.xcorrFrequency = 1000;  % unit in Hz
s.markerSize = 500;
s.pawNames = {'LH', 'LF', 'RF', 'RH'};
s.pawColors = hsv(4);
s.unit_id = [];

% reassign settings contained in opts
if exist('opts', 'var'); for i = 1:2:length(opts); s.(opts{i}) = opts{i+1}; end; end

% check if the folder for saving figs exists
if ~exist(s.folder, 'dir')
    mkdir(s.folder)
end

% loading files
disp('loading runAnalyzed.mat...')
sessionFolder = fullfile(getenv('OBSDATADIR'), 'sessions', session);
load(fullfile(sessionFolder, 'runAnalyzed.mat'));

disp('loading kinData.mat...');
if exist(fullfile(sessionFolder, 'kinData.mat'))
    load(fullfile(sessionFolder, 'kinData.mat'), 'kinData', 'stanceBins')
else
    [kinData, stanceBins] = getKinematicData(session);
end

disp('loading neuralData.mat...');
load(fullfile(sessionFolder, 'neuralData.mat'));


for unit_ind = 1:length(unit_ids)
    fprintf('%s: plotting cell %i/%i\n', session, unit_ind, length(unit_ids))
    unit_spkRates = spkRates(unit_ind, :);
    
    if ~isempty(s.unit_id)
        unit_ind = length(unit_ids);
        unit_spkRates = spkRates(find(unit_ids == unit_id), :);
        unit_id = s.unit_id;
    else
        unit_spkRates = spkRates(unit_ind, :);
        unit_id = unit_ids(unit_ind);
    end
    
    
    % get rid of nans in neural firing rates
    startInd = find(~isnan(unit_spkRates), 1, 'first');
    endInd = find(~isnan(unit_spkRates), 1, 'last');
    unit_spkRates = unit_spkRates(startInd:endInd);
    neuralTimes = timeStamps(startInd:endInd);
    
    % restricting wheel times and vel by start and end of neural times
    vel = getVelocity(wheelPositions, .02, 1000); % should have the same length as wheelTimes
    ind(1) = knnsearch(wheelTimes', neuralTimes(1));
    ind(2) = knnsearch(wheelTimes', neuralTimes(end));
    velTimeStamps = wheelTimes(ind(1) : ind(2));
    vel = vel(ind(1) : ind(2));
    
    
    % loading paw and tail positions
    disp('getting paw and tail positions...');
    locationsTable = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw.csv')); % get raw tracking data
    [locations, features] = fixTracking(locationsTable, frameTimeStamps, pixelsPerM);
    
    % restricting frames by neural times
    startTimeInd = knnsearch(frameTimeStamps, neuralTimes(1));
    endTimeInd = knnsearch(frameTimeStamps, neuralTimes(end));
    selectedFrameTimeStamps = frameTimeStamps(startTimeInd : endTimeInd);
    
    % get velocity for each frame
    disp('getting frame vels..');
    if size(velTimeStamps, 1) == 1
        velTimeStamps = velTimeStamps';
    end
    frameVel = nan(length(selectedFrameTimeStamps), 1);
    
    % find corresponding velocity for every frame times
    temp = selectedFrameTimeStamps(:, 1);
    [~,ind] = histc(temp,[-inf;(velTimeStamps(1:end-1)+velTimeStamps(2:end))/2;inf]);
    selectedFrameTimeStamps(:, 2) = velTimeStamps(ind(:, 1));
    frameVel = vel(ind(:, 1))';
    
    [runningBins, runningStartInds, runningEndInds] = isRunning(frameVel);
    
    % processing paws x y z positions
    botPawInds = find(contains(features, 'paw') & contains(features, '_bot'));
    topPawInds = find(contains(features, 'paw') & contains(features, '_top'));
    
    locationsPaws = nan(size(locations,1), 3, 4);  % frameNum X xyz X pawNum
    locationsPaws(:,1:2,:) = locations(:,:,botPawInds);
    locationsPaws(:,3,:) = locations(:,2,topPawInds);
    locationsPawsPixels = locationsPaws;  % copy version in pixel coordinates before making subsequent transformations
    locationsPaws(:,2,:) = locationsPaws(:,2,:) - nosePos(2); % subtract midline from all y values
    locationsPaws(:,3,:) = (wheelCenter(2)-wheelRadius) - locationsPaws(:,3,:); % flip z and set s.t. top of wheel is zero
    
    locationsPaws = locationsPaws/pixelsPerM;
    pawsXLocations = squeeze(locationsPaws(startTimeInd:endTimeInd, 1, :));
    pawsYLocations = squeeze(locationsPaws(startTimeInd:endTimeInd, 2, :));
    pawsZLocations = squeeze(locationsPaws(startTimeInd:endTimeInd, 3, :));
    
    % processing tail x y z positions
    botTailBins = contains(features, 'tail') & contains(features, '_bot');
    topTailBins = contains(features, 'tail') & contains(features, '_top');
    locationsTail = nan(size(locations,1), 3, 2);  % frameNum X xyz X tailbase/mid
    locationsTail(:,1:2,:) = locations(:,:,botTailBins);
    locationsTail(:,3,:) = locations(:,2,topTailBins);
    locationsTail(:,2,:) = locationsTail(:,2,:) - nosePos(2); % subtract midline from all y values
    locationsTail(:,3,:) = (wheelCenter(2)-wheelRadius) - locationsTail(:,3,:); % flip z and set s.t. top of wheel is zero
    
    locationsTail = locationsTail/pixelsPerM;
    locationsTail = locationsTail(startTimeInd:endTimeInd, :, :);
    
    % -----------------------------------------------------------------------%
    % Figures paw position
    figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
    plotInd = 1;
    for i = 1:4
        plotInd = plotContinuousData(pawsZLocations(:, i), unit_spkRates, plotInd, logical(runningBins), selectedFrameTimeStamps, neuralTimes, ...
            {'dataName', [s.pawNames{i} ' paw z positions']});
    end
    fileName = fullfile(s.folder, ['unit' num2str(unit_id) '_paws_Zpositions']);
    saveas(gcf, [fileName '.png'])
    
    
    figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
    plotInd = 1;
    for i = 1:4
        plotInd = plotContinuousData(pawsXLocations(:, i), unit_spkRates, plotInd, logical(runningBins), selectedFrameTimeStamps, neuralTimes, ...
            {'dataName', [s.pawNames{i} ' paw x positions']});
    end
    fileName = fullfile(s.folder, ['unit' num2str(unit_id) '_paws_Xpositions']);
    saveas(gcf, [fileName '.png'])
    
    
    figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
    plotInd = 1;
    for i = 1:4
        plotInd = plotContinuousData(pawsYLocations(:, i), unit_spkRates, plotInd, logical(runningBins), selectedFrameTimeStamps, neuralTimes, ...
            {'dataName', [s.pawNames{i} ' paw y positions']});
    end
    fileName = fullfile(s.folder, ['unit' num2str(unit_id) '_paws_Ypositions']);
    saveas(gcf, [fileName '.png'])
    
    % ----------------------------------------------------------------------- %
    % Figure velocity in x y z axis
    
    figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
    plotInd = 1;
    frameFreq = 250; % Hz
    xVel = nan(length(selectedFrameTimeStamps)-1, 4);
    for i = 1:4
        xVel(:, i) = abs(diff(pawsXLocations(:, i))) / (1/frameFreq);
        plotInd = plotContinuousData(xVel(:, i), unit_spkRates, plotInd, logical(runningBins), selectedFrameTimeStamps(2:end), neuralTimes, ...
            {'dataName', [s.pawNames{i} ' paw x velocity (m/s)']});
    end
    fileName = fullfile(s.folder, ['unit' num2str(unit_id) '_paws_Xvelocity']);
    saveas(gcf, [fileName '.png'])
    
    figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
    plotInd = 1;
    yVel = nan(length(selectedFrameTimeStamps)-1, 4);
    for i = 1:4
        yVel(:, i) = abs(diff(pawsYLocations(:, i))) / (1/frameFreq);
        plotInd = plotContinuousData(yVel(:, i), unit_spkRates, plotInd, logical(runningBins), selectedFrameTimeStamps(2:end), neuralTimes, ...
            {'dataName', [s.pawNames{i} ' paw y velocity (m/s)']});
    end
    fileName = fullfile(s.folder, ['unit' num2str(unit_id) '_paws_Yvelocity']);
    saveas(gcf, [fileName '.png'])
    
    
    figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
    plotInd = 1;
    zVel = nan(length(selectedFrameTimeStamps)-1, 4);
    for i = 1:4
        zVel(:, i) = abs(diff(pawsZLocations(:, i))) / (1/frameFreq);
        plotInd = plotContinuousData(zVel(:, i), unit_spkRates, plotInd, logical(runningBins), selectedFrameTimeStamps(2:end), neuralTimes, ...
            {'dataName', [s.pawNames{i} ' paw z velocity (m/s)']});
    end
    fileName = fullfile(s.folder, ['unit' num2str(unit_id) '_paws_Zvelocity']);
    saveas(gcf, [fileName '.png'])
    
    % ----------------------------------------------------------------------- %
    % Figure acceleration in x y z axis
    
    figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
    plotInd = 1;
    xAcc = nan(length(selectedFrameTimeStamps)-2, 4);
    
    for i = 1:4
        xAcc(:, i) = abs(diff(xVel(:, i))) / (1/frameFreq);
        plotInd = plotContinuousData(xAcc(:, i), unit_spkRates, plotInd, logical(runningBins), selectedFrameTimeStamps(3:end), neuralTimes, ...
            {'dataName', [s.pawNames{i} ' paw x Acceleration (m/s^2)']});
    end
    fileName = fullfile(s.folder, ['unit' num2str(unit_id) '_paws_Xacceleration']);
    saveas(gcf, [fileName '.png'])
    
    
    figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
    plotInd = 1;
    yAcc = nan(length(selectedFrameTimeStamps)-2, 4);
    
    for i = 1:4
        yAcc(:, i) = abs(diff(yVel(:, i))) / (1/frameFreq);
        plotInd = plotContinuousData(yAcc(:, i), unit_spkRates, plotInd, logical(runningBins), selectedFrameTimeStamps(3:end), neuralTimes, ...
            {'dataName', [s.pawNames{i} ' paw y Acceleration (m/s^2)']});
    end
    fileName = fullfile(s.folder, ['unit' num2str(unit_id) '_paws_Yacceleration']);
    saveas(gcf, [fileName '.png'])
    
    figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
    plotInd = 1;
    zAcc = nan(length(selectedFrameTimeStamps)-2, 4);
    
    for i = 1:4
        zAcc(:, i) = abs(diff(zVel(:, i))) / (1/frameFreq);
        plotInd = plotContinuousData(zAcc(:, i), unit_spkRates, plotInd, logical(runningBins), selectedFrameTimeStamps(3:end), neuralTimes, ...
            {'dataName', [s.pawNames{i} ' paw z Acceleration (m/s^2)']});
    end
    fileName = fullfile(s.folder, ['unit' num2str(unit_id) '_paws_Zacceleration']);
    saveas(gcf, [fileName '.png'])
    
    % ----------------------------------------------------------------------- %
    % Figures velocity + body angles + tail angles
    
    figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
    plotInd = 1;
    % velocity
    plotInd = plotContinuousData(frameVel, unit_spkRates, plotInd, logical(runningBins), selectedFrameTimeStamps, neuralTimes, ...
        {'dataName', 'velocity (m/s)'});
    % bodyAngle
    plotInd = plotContinuousData(bodyAngles, unit_spkRates, plotInd, logical(runningBins), selectedFrameTimeStamps, neuralTimes, ...
        {'dataName', 'bodyAngles'});
    
    % tailAngle
    tailBase = squeeze(locationsTail(:, :, 1));
    tailMid = squeeze(locationsTail(:, :, 2));
    
    tailBase = fillmissing(tailBase, 'spline');
    tailMid = fillmissing(tailMid, 'spline');
    
    tailAngles = (tailBase(:, 2) - tailMid(:, 2)) ./ (tailBase(:, 1) - tailMid(:, 1));
    
    plotInd = plotContinuousData(tailAngles, unit_spkRates, plotInd, logical(runningBins), selectedFrameTimeStamps, neuralTimes, ...
        {'dataName', 'tailAngles'});
    
    
    % body z angles
    noseZPositions = locations(startTimeInd:endTimeInd, :, 7);
    noseZPositions(:, 2) = (wheelCenter(2)-wheelRadius) - noseZPositions(:, 2);
    noseZPositions = noseZPositions / pixelsPerM;
    
    tailBaseTop = locations(startTimeInd:endTimeInd, :, 5);
    tailBaseTop(:, 2) = (wheelCenter(2)-wheelRadius) - tailBaseTop(:, 2);
    tailBaseTop = tailBaseTop / pixelsPerM;
    tailBaseTop = fillmissing(tailBaseTop, 'spline');
    
    bodyZAngles = (noseZPositions(:, 2) - tailBaseTop(:, 2)) ./ (noseZPositions(:, 1) - tailBaseTop(:, 1));
    
    plotInd = plotContinuousData(bodyZAngles, unit_spkRates, plotInd, logical(runningBins), selectedFrameTimeStamps, neuralTimes, ...
        {'dataName', 'bodyZAngles'});
    
    fileName = fullfile(s.folder, ['unit' num2str(unit_id) '_vel_bodyAngles_tailAngles']);
    saveas(gcf, [fileName '.png']);
    
    disp('Unit Done!');
    close all
end

disp('Session DONE!');


end
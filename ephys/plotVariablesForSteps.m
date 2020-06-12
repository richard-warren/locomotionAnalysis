function plotVariablesForSteps(session, opts)


% settings
s.folder = fullfile(getenv('OBSDATADIR'), 'figures', 'ephys', 'VariablesForSteps', session);  % folder in which the plots will be saved
s.stepPercentiles = [40 60]; % only include steps with durations in between these percentile limits
s.unit_id = 17;
s.pawNames = {'LH', 'LF', 'RF', 'RH'};
s.rows = 4;
s.cols = 3;
s.fontSize = 12;
s.lineWidth = 3;

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

% if users ask to plot for all the units in one session
if strcmp(s.unit_id, 'all')
    s.unit_id = unit_ids;
end
    

for unitInd = 1:length(s.unit_id)
    fprintf('%s: plotting cell %i/%i\n', session, unitInd, length(s.unit_id))
    
    unit = s.unit_id(unitInd);
    unit_spkRates = spkRates(find(unit_ids == unit), :);
    
    % get rid of nans in neural firing rates
    disp('reformatting neural data...');
    neuralStartTimeInd = find(~isnan(unit_spkRates), 1, 'first');
    neuralEndTimeInd = find(~isnan(unit_spkRates), 1, 'last');
    unit_spkRates = unit_spkRates(neuralStartTimeInd : neuralEndTimeInd);
    neuralTimes = timeStamps(neuralStartTimeInd : neuralEndTimeInd);
    
    % restricting wheel times and vel by start and end of neural times
    disp('getting vel...');
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
    stanceBins = stanceBins(startTimeInd:endTimeInd, :);

    
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
    
%     % processing tail x y z positions
%     botTailBins = contains(features, 'tail') & contains(features, '_bot');
%     topTailBins = contains(features, 'tail') & contains(features, '_top');
%     locationsTail = nan(size(locations,1), 3, 2);  % frameNum X xyz X tailbase/mid
%     locationsTail(:,1:2,:) = locations(:,:,botTailBins);
%     locationsTail(:,3,:) = locations(:,2,topTailBins);
%     locationsTail(:,2,:) = locationsTail(:,2,:) - nosePos(2); % subtract midline from all y values
%     locationsTail(:,3,:) = (wheelCenter(2)-wheelRadius) - locationsTail(:,3,:); % flip z and set s.t. top of wheel is zero
%     
%     locationsTail = locationsTail/pixelsPerM;
%     locationsTail = locationsTail(startTimeInd:endTimeInd, :, :);
    
    phase = struct();
    neuralA = struct();
    positions = struct();
    zpositions = struct();
    for i = 1:4
        
        % get swing start and end times
        swingStartInds = find(diff(~stanceBins(:, i))==1);
        swingStartTimes = selectedFrameTimeStamps(swingStartInds(1:end-1));
        stanceEndTimes = selectedFrameTimeStamps(swingStartInds(2:end)-1);
        sessionEvents = cat(2, swingStartTimes, stanceEndTimes);
        sessionEvents = sessionEvents(~isnan(sum(sessionEvents,2)),:); % remove nan entries
        
        % only take steps in middle of duration distribution
        durations = diff(sessionEvents,1,2);
        durationLimits = prctile(durations, s.stepPercentiles);
        inds = find(durations>durationLimits(1) & durations<durationLimits(2));
        sessionEvents = sessionEvents(durations>durationLimits(1) & durations<durationLimits(2), :);        
        
        stanceEndInds = swingStartInds(inds+1)-1;
        swingStartInds = swingStartInds(inds);        
       
        for j = 1:length(sessionEvents)
            x = selectedFrameTimeStamps(swingStartInds(j):stanceEndInds(j));
            positionsTemp = pawsXLocations(swingStartInds(j):stanceEndInds(j), i);
            yPositioinsTemp = pawsYLocations(swingStartInds(j):stanceEndInds(j), i);
            zPositionsTemp = pawsZLocations(swingStartInds(j):stanceEndInds(j), i);
            positionsTemp = positionsTemp - mean(positionsTemp);
            yPositioinsTemp = yPositioinsTemp - mean(yPositioinsTemp);
            zPositionsTemp = zPositionsTemp - mean(zPositionsTemp);
            positions.(s.pawNames{i})(j, :) = interp1(1:length(positionsTemp), positionsTemp, linspace(1, length(positionsTemp), 100));
            zpositions.(s.pawNames{i})(j, :) = interp1(1:length(zPositionsTemp), zPositionsTemp, linspace(1, length(zPositionsTemp), 100));
            ypositions.(s.pawNames{i})(j, :) = interp1(1:length(yPositioinsTemp), yPositioinsTemp, linspace(1, length(yPositioinsTemp), 100));
            % positions = lowpass(positions, 5, 250);
            hPositions = hilbert(positionsTemp);
            y = angle(hPositions) + pi;
            phase.(s.pawNames{i})(j, :) = interp1(x, y, linspace(x(1), x(end), 100));
            
            bins = neuralTimes > sessionEvents(j, 1) & neuralTimes < sessionEvents(j, 2);
            nx = neuralTimes(bins);
            ny = unit_spkRates(bins);
            neuralA.(s.pawNames{i})(j, :) = interp1(nx, ny, linspace(nx(1), nx(end), 100));
        end
        
        % plot avg traces for each variable
        disp(['plotting avg traces...paw ' num2str(i)]);
        figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
        reformatFun = @(x) nanmean((x - nanmean(x, 2))./max(abs(x - nanmean(x, 2)), [], 2));
        plot(reformatFun(phase.(s.pawNames{i})), 'LineWidth', 1);
        plot(reformatFun(neuralA.(s.pawNames{i})), 'LineWidth', 1);
        plot(reformatFun(zpositions.(s.pawNames{i})), 'LineWidth', 1);
        plot(reformatFun(positions.(s.pawNames{i})), 'LineWidth', 1);
        plot(reformatFun(ypositions.(s.pawNames{i})), 'LineWidth', 1);
        xlabel('swing start -> stance end');
        legend([s.pawNames{i} ' limb phase'], 'neural FR', [s.pawNames{i} ' limb z position'], [s.pawNames{i} ' limb x position'], [s.pawNames{i} ' limb y position']);
     
        fileName = fullfile(s.folder, ['unit' num2str(unit) '_variablesStepsAvg_' s.pawNames{i}]);
        saveas(gcf, [fileName '.png'])    
        
    end
    
    close all
    
 
    
    disp('plotting cross corr..');
    figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
    plotInd = 1;

    for paw = 1:4 
        c = struct();
        for i = 1:length(neuralA.(s.pawNames{paw}))
            x = positions.(s.pawNames{paw})(i, :);
            y = ypositions.(s.pawNames{paw})(i, :);
            z = zpositions.(s.pawNames{paw})(i, :);
            neuA = neuralA.(s.pawNames{paw})(i, :);
            [c.x(i, :), lags] = xcorr(x - mean(x), neuA - mean(neuA), 30, 'normalized');
            [c.y(i, :), lags] = xcorr(y - mean(y), neuA - mean(neuA), 30, 'normalized');
            [c.z(i, :), lags] = xcorr(z - mean(z), neuA - mean(neuA), 30, 'normalized');        
        end
    
        
        % x positions
        subplot(s.rows, s.cols, plotInd); 
        xrange = linspace(-30, 30, size(c.x, 2));
        plot(xrange, nanmean(c.x, 1), '-', 'LineWidth', s.lineWidth, 'Color', [0 0.4470 0.7410]); hold on;
        
        minind = find(nanmean(c.x, 1) == min(nanmean(c.x, 1)));
        xlimit = get(gca,'xlim');
        ylimit = get(gca,'ylim');
        text(xrange(minind), min(nanmean(c.x, 1)), ['min corr = ' num2str(min(nanmean(c.x, 1)))], 'FontSize', s.fontSize);
        text(xrange(minind), min(nanmean(c.x, 1)) + 0.1*(ylimit(2) - ylimit(1)), ['min.x = ' num2str(xrange(minind))], 'FontSize', s.fontSize);
        line([xrange(minind), xrange(minind)], [ylimit(1), ylimit(2)], 'Color', [1 0.77 0.16], 'LineWidth', s.lineWidth);
        
        maxind = find(nanmean(c.x, 1) == max(nanmean(c.x, 1)));
        text(xrange(maxind), max(nanmean(c.x, 1)), ['max corr = ' num2str(max(nanmean(c.x, 1)))], 'FontSize', s.fontSize);
        text(xrange(maxind), max(nanmean(c.x, 1)) - 0.1*(ylimit(2) - ylimit(1)), ['max.x = ' num2str(xrange(maxind))], 'FontSize', s.fontSize);
        line([xrange(maxind), xrange(maxind)], [ylimit(1), ylimit(2)], 'Color', [1 0.42 0.54], 'LineWidth', s.lineWidth);
        
      
        ylabel('corr coef');
        xlabel('% of one step cycle')
        title(['corss corr: ', s.pawNames{paw}, ' paw and x positions'], 'FontSize', s.fontSize);
        plotInd = plotInd + 1;
        
        % y positions
        subplot(s.rows, s.cols, plotInd); 
        xrange = linspace(-30, 30, size(c.y, 2));
        plot(xrange, nanmean(c.y, 1), '-', 'LineWidth', s.lineWidth, 'Color', [0 0.4470 0.7410]); hold on;
        
        minind = find(nanmean(c.y, 1) == min(nanmean(c.y, 1)));
        xlimit = get(gca,'xlim');
        ylimit = get(gca,'ylim');
        text(xrange(minind), min(nanmean(c.y, 1)), ['min corr = ' num2str(min(nanmean(c.y, 1)))], 'FontSize', s.fontSize);
        text(xrange(minind), min(nanmean(c.y, 1)) + 0.1*(ylimit(2) - ylimit(1)), ['min.x = ' num2str(xrange(minind))], 'FontSize', s.fontSize);
        line([xrange(minind), xrange(minind)], [ylimit(1), ylimit(2)], 'Color', [1 0.77 0.16], 'LineWidth', s.lineWidth);
        
        maxind = find(nanmean(c.y, 1) == max(nanmean(c.y, 1)));
        text(xrange(maxind), max(nanmean(c.y, 1)), ['max corr = ' num2str(max(nanmean(c.y, 1)))], 'FontSize', s.fontSize);
        text(xrange(maxind), max(nanmean(c.y, 1)) - 0.1*(ylimit(2) - ylimit(1)), ['max.x = ' num2str(xrange(maxind))], 'FontSize', s.fontSize);
        line([xrange(maxind), xrange(maxind)], [ylimit(1), ylimit(2)], 'Color', [1 0.42 0.54], 'LineWidth', s.lineWidth);
        
        
        ylabel('corr coef');
        xlabel('% of one step cycle')
        title(['corss corr: ', s.pawNames{paw}, ' paw and y positions'], 'FontSize', s.fontSize);
        plotInd = plotInd + 1;
        
        
        % z positions
        subplot(s.rows, s.cols, plotInd); 
        xrange = linspace(-30, 30, size(c.z, 2));
        plot(xrange, nanmean(c.z, 1), '-', 'LineWidth', s.lineWidth, 'Color', [0 0.4470 0.7410]); hold on;
        
        minind = find(nanmean(c.z, 1) == min(nanmean(c.z, 1)));
        xlimit = get(gca,'xlim');
        ylimit = get(gca,'ylim');
        text(xrange(minind), min(nanmean(c.z, 1)), ['min corr = ' num2str(min(nanmean(c.z, 1)))], 'FontSize', s.fontSize);
        text(xrange(minind), min(nanmean(c.z, 1)) + 0.1*(ylimit(2) - ylimit(1)), ['min.x = ' num2str(xrange(minind))], 'FontSize', s.fontSize);
        line([xrange(minind), xrange(minind)], [ylimit(1), ylimit(2)], 'Color', [1 0.77 0.16], 'LineWidth', s.lineWidth);
        
        maxind = find(nanmean(c.z, 1) == max(nanmean(c.z, 1)));
        text(xrange(maxind), max(nanmean(c.z, 1)), ['max corr = ' num2str(max(nanmean(c.z, 1)))], 'FontSize', s.fontSize);
        text(xrange(maxind), max(nanmean(c.z, 1)) - 0.1*(ylimit(2) - ylimit(1)), ['max.x = ' num2str(xrange(maxind))], 'FontSize', s.fontSize);
        line([xrange(maxind), xrange(maxind)], [ylimit(1), ylimit(2)], 'Color', [1 0.42 0.54], 'LineWidth', s.lineWidth);
        
      
        ylabel('corr coef');
        xlabel('% of one step cycle')
        title(['corss corr: ', s.pawNames{paw}, ' paw and z positions'], 'FontSize', s.fontSize);
        plotInd = plotInd + 1; 
        
    end
    
    fileName = fullfile(s.folder, ['unit' num2str(unit) '_variablesCrossCorr']);
    saveas(gcf, [fileName '.png']) 
    
    disp('Unit Done!')   

end
  
disp('Session Done!')


end
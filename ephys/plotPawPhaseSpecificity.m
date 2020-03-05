function plotPawPhaseSpecificity(session, unit_id, paw_id, opts)

% settings
s.folder = fullfile(getenv('OBSDATADIR'), 'figures', 'ephys', 'pawPhaseSpecificity');  % folder in which the plots will be saved
s.phaseMode = 'stance'; % indicate the plot should focus on swing start time points or stance start time points
s.stepPercentiles = [30 70]; % only include steps with durations in between these percentile limits
s.pawNames = {'LH', 'LF', 'RF', 'RH'};
s.pawColors = hsv(4);
s.mainColor = [.8 .4 1];
s.markerSize = 20;

% reassign settings contained in opts
if exist('opts', 'var'); for i = 1:2:length(opts); s.(opts{i}) = opts{i+1}; end; end

% initializations
sessionFolder = fullfile(getenv('OBSDATADIR'), 'sessions', session);
disp('loading runAnalyzed.mat...');
load(fullfile(sessionFolder, 'runAnalyzed.mat'), ...
        'obsOnTimes', 'obsOffTimes',  'wiskContactFrames', 'frameTimeStamps', ...
        'frameTimeStampsWisk', 'rewardTimes', 'isLightOn', 'touches', 'touchesPerPaw', 'touchClassNames', 'wheelPositions', 'wheelTimes');
vel = getVelocity(wheelPositions, .02, 1000); % should have the same length as wheelTimes
velTimeStamps = wheelTimes;
clear wheelTimes
clear wheelPositions

% load kinData if it exists // otherwise compute kinData    
disp('loading kinData...');
if exist(fullfile(sessionFolder, 'kinData.mat'))
    load(fullfile(sessionFolder, 'kinData.mat'), 'kinData', 'stanceBins')
else
    [kinData, stanceBins] = getKinematicData(session);
end

% load neuralData
disp('loading neuralData.mat...');
load(fullfile(sessionFolder, 'neuralData.mat'));
spkTimes = spkTimes';
unit_spkTimes = spkTimes{find(unit_ids == unit_id), 1};  % an N-by-1 matrix

% load cellData
disp('loading cellData.csv...');
cellData = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'cellData.csv'));
unit_idInd = find(cellData.unit_id == unit_id);
unit_startTime = polyval(openEphysToSpikeMapping, cellData.timeStart(unit_idInd)*60);
unit_endTime = polyval(openEphysToSpikeMapping, cellData.timeEnd(unit_idInd)*60);

targetPaw = paw_id;

switch paw_id
    case 1
        diagPaw = 3;
    case 2
        diagPaw = 4;
    case 3
        diagPaw = 1;
    case 4
        diagPaw = 2;
end


% get swing start and end times for targer paw
swingStartInds = find(diff(~stanceBins(:,targetPaw))==1);
swingStartTimes = frameTimeStamps(swingStartInds(1:end-1));
stanceStartTimes = frameTimeStamps(swingStartInds(2:end)-1);
swingStanceTimes = cat(2, swingStartTimes, stanceStartTimes);
swingStanceTimes = swingStanceTimes(~isnan(sum(swingStanceTimes,2)),:); % remove nan entries

% only take steps in middle of duration distribution
durations = swingStanceTimes(:, 2) - swingStanceTimes(:, 1);
durationLimits = prctile(durations, s.stepPercentiles);
validInds = find(durations>durationLimits(1) & durations<durationLimits(2));
validSwingStartTimes = swingStartTimes(validInds, :);
validStanceStartTimes = stanceStartTimes(validInds, :);
validSwingStanceTimes = swingStanceTimes(validInds, :);

% bin steps by unit spk time
validInds = find(validSwingStartTimes >= unit_startTime & validStanceStartTimes <= unit_endTime);
validSwingStartTimes = swingStartTimes(validInds, :);
validStanceStartTimes = stanceStartTimes(validInds, :);
validSwingStanceTimes = swingStanceTimes(validInds, :);

% bin steps by velocity
validInds = [];
for i = 1:length(validSwingStartTimes)-1
    avgVel = mean(vel(velTimeStamps >= validSwingStartTimes(i) & velTimeStamps <= validSwingStartTimes(i+1)));
    if avgVel >= 0.3
        validInds = [validInds; i];
    end
end

validSwingStartTimes = swingStartTimes(validInds, :);
validStanceStartTimes = stanceStartTimes(validInds, :);
validSwingStanceTimes = swingStanceTimes(validInds, :);





% get swing start and end times for the diagnal paw, which is in phase with
% the target paw
swingStartInds_d = find(diff(~stanceBins(:,diagPaw))==1);
swingStartTimes_d = frameTimeStamps(swingStartInds_d(1:end-1));
stanceStartTimes_d = frameTimeStamps(swingStartInds_d(2:end)-1);
swingStanceTimes_d = cat(2, swingStartTimes_d, stanceStartTimes_d);
swingStanceTimes_d = swingStanceTimes_d(~isnan(sum(swingStanceTimes_d,2)),:); % remove nan entries

validInds_d = knnsearch(swingStanceTimes_d(:, 1), validSwingStanceTimes(:, 1));
validSwingStanceTimes_d = swingStanceTimes_d(validInds_d, :);


switch s.phaseMode
    case 'stance'
        % get the phase offset b/w target limb and its diagnal limb
        phaseOffset = validSwingStanceTimes(:, 2) - validSwingStanceTimes_d(:, 2);
        phaseOffset = phaseOffset(1:end-1);
        [stepPhaseOffset, stepID] = sort(phaseOffset);
        temp = [validSwingStartTimes(1:end-1), validStanceStartTimes(1:end-1), validSwingStartTimes(2:end)];
        temp = temp - validStanceStartTimes(1:end-1);
        
        % get the spike times for every valid swing stance cycle
        % and PLOT!!!
        
        fprintf('%s: plotting unit %i', session, unit_id)
        figure('name', sprintf('%s - unit %i', session, unit_id), ...
        'color', 'white', 'MenuBar', 'none', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
        
        for i = 1:length(validSwingStartTimes)-1
            stepSpikeTimes = unit_spkTimes(unit_spkTimes > validSwingStanceTimes(i, 1) & unit_spkTimes < validSwingStanceTimes(i+1, 1));
            xValue = stepSpikeTimes - validSwingStanceTimes(i, 2);
            y = stepPhaseOffset(find(stepID == i)); % unit:second
            yValue = repmat(y, size(xValue, 1), 1);
            
            scatter(xValue, yValue, s.markerSize, '.k');
       
%             xlabel('stance start');
%             ylabel('diagonal paw phase offset');
            
            
        end

end



end
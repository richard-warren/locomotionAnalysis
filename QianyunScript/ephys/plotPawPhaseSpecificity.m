function plotPawPhaseSpecificity(session, unit_id, paw_id, opts)

% settings
s.folder = fullfile(getenv('OBSDATADIR'), 'figures', 'ephys', 'pawPhaseSpecificity');  % folder in which the plots will be saved
s.phaseMode = 'stance'; % indicate the plot should focus on swing start time points or stance start time points
s.stepPercentiles = [30 70]; % only include steps with durations in between these percentile limits
s.plotStepWindow = 0.2; % unit in second. 
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
swingStartTimes = frameTimeStamps(swingStartInds);
validBins = find(swingStartTimes > unit_startTime & swingStartTimes < unit_endTime);
validSwingTimes = swingStartTimes(validBins);



stanceStartInds = find(diff(stanceBins(:,targetPaw))== 1);
stanceStartTimes = frameTimeStamps(stanceStartInds);
validBins = find(stanceStartTimes > unit_startTime & stanceStartTimes < unit_endTime);
validStanceTimes = stanceStartTimes(validBins);


if validSwingTimes(1) < validStanceTimes(1)   
    if length(validSwingTimes) ~= length(validStanceTimes)
        validSwingTimes(end) = [];
    end    
    
else 
    validStanceTimes(1) = [];
    
    if length(validSwingTimes) ~= length(validStanceTimes)
        validSwingTimes(end) = [];
    end
end  
     

validSwingStanceTimes = cat(2, validSwingTimes, validStanceTimes);
validSwingStanceTimes = validSwingStanceTimes(~isnan(sum(validSwingStanceTimes,2)),:); % remove nan entries





% bin steps by velocity
validInds = [];
for i = 1:length(validSwingTimes)-1
    avgVel = mean(vel(velTimeStamps >= validSwingTimes(i) & velTimeStamps <= validSwingTimes(i+1)));
    if avgVel >= 0.3
        validInds = [validInds; i];
    end
end

validSwingStanceTimes = validSwingStanceTimes(validInds, :);




% only take steps in middle of duration distribution
durations = validSwingStanceTimes(:, 2) - validSwingStanceTimes(:, 1);
durationLimits = prctile(durations, s.stepPercentiles);
validInds = find(durations>durationLimits(1) & durations<durationLimits(2));
validSwingStanceTimes = validSwingStanceTimes(validInds, :);


% only take durations that are less than 0.5s
durations = validSwingStanceTimes(:, 2) - validSwingStanceTimes(:, 1);
validInds = find(durations < 0.5);
validSwingStanceTimes = validSwingStanceTimes(validInds, :);


 % get swing start and end times for the diagnal paw, which is in phase with the target paw

swingStartInds_d = find(diff(~stanceBins(:,diagPaw))==1);
swingStartTimes_d = frameTimeStamps(swingStartInds_d);
validBins = find(swingStartTimes_d > unit_startTime & swingStartTimes_d < unit_endTime);
validSwingTimes_d = swingStartTimes_d(validBins);

stanceStartInds_d = find(diff(stanceBins(:,diagPaw))== 1);
stanceStartTimes_d = frameTimeStamps(stanceStartInds_d);
validBins = find(stanceStartTimes_d > unit_startTime & stanceStartTimes_d < unit_endTime);
validStanceTimes_d = stanceStartTimes_d(validBins);

if validSwingTimes_d(1) < validStanceTimes_d(1)   
    if length(validSwingTimes_d) ~= length(validStanceTimes_d)
        validSwingTimes_d(end) = [];
    end    
    
else 
    validStanceTimes_d(1) = [];
    
    if length(validSwingTimes_d) ~= length(validStanceTimes_d)
        validSwingTimes_d(end) = [];
    end
end 

validSwingStanceTimes_d = cat(2, validSwingTimes_d, validStanceTimes_d);
validSwingStanceTimes_d = validSwingStanceTimes_d(~isnan(sum(validSwingStanceTimes_d,2)),:); % remove nan entries





% Plot the raster!

switch s.phaseMode
    case 'stance'       
        
        % get matched steps for the diagnal paw
        validInds_d = knnsearch(validSwingStanceTimes_d(:, 2), validSwingStanceTimes(:, 2));
        validSwingStanceTimes_d = validSwingStanceTimes_d(validInds_d, :);
        
        
        % exclude steps that diagnal paw duration is bigger than 0.5s
        durations_d = validSwingStanceTimes_d(:, 2) - validSwingStanceTimes_d(:, 1);
        validInds = find(durations_d < 0.5);
        validSwingStanceTimes = validSwingStanceTimes(validInds, :);
        validSwingStanceTimes_d = validSwingStanceTimes_d(validInds, :);
               
        
        % plot #1: focus on the target paw
        fprintf('%s: plotting unit %i', session, unit_id)
        figure('name', sprintf('%s - unit %i', session, unit_id), ...
            'color', 'white', 'MenuBar', 'none', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
               
        % get the phase offset b/w target limb and its diagnal limb
        phaseOffset = validSwingStanceTimes(:, 2) - validSwingStanceTimes_d(:, 2);
        % get times for plotting spikes        
        times = [validSwingStanceTimes(:, 2) - s.plotStepWindow, validSwingStanceTimes(:, 2) + s.plotStepWindow];
        
        
        for i = 1:length(validSwingStanceTimes)
            stepSpikeTimes = unit_spkTimes(unit_spkTimes > times(i, 1) & unit_spkTimes < times(i, 2));
            xValue = stepSpikeTimes - validSwingStanceTimes(i, 2);
            y = phaseOffset(i); % unit:second
            yValue = repmat(y, size(xValue, 1), 1);
            
            scatter(xValue, yValue, s.markerSize, '.k');
            
%             s.XLabel = 'target paw stance start time';
%             s.YLabel = 'diagonal paw stance time offset';
                      
        end
        % save
        fileName = fullfile(s.folder, strcat(session, '_unit', num2str(unit_id), '_paw_', s.pawNames{paw_id}));
        % savefig(fileName)      
        saveas(gcf, strcat(fileName, '.png'))
        disp('file saved');
        
        
        
         % plot #2: focus on the diagnal paw
        fprintf('%s: plotting unit %i', session, unit_id)
        figure('name', sprintf('%s - unit %i - diagnal paw', session, unit_id), ...
            'color', 'white', 'MenuBar', 'none', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
               
        % get the phase offset b/w target limb and its diagnal limb
        phaseOffset = validSwingStanceTimes_d(:, 2) -  validSwingStanceTimes(:, 2);
        % get times for plotting spikes        
        times = [validSwingStanceTimes_d(:, 2) - s.plotStepWindow, validSwingStanceTimes_d(:, 2) + s.plotStepWindow];
        
        
        for i = 1:length(validSwingStanceTimes)
            stepSpikeTimes = unit_spkTimes(unit_spkTimes > times(i, 1) & unit_spkTimes < times(i, 2));
            xValue = stepSpikeTimes - validSwingStanceTimes_d(i, 2);
            y = phaseOffset(i); % unit:second
            yValue = repmat(y, size(xValue, 1), 1);
            
            scatter(xValue, yValue, s.markerSize, '.k');
            
%             s.XLabel = 'diagnal paw stance start time';
%             s.YLabel = 'target paw stance time offset';
                      
        end
        % save
        fileName = fullfile(s.folder, strcat(session, '_unit', num2str(unit_id), '_paw_', s.pawNames{diagPaw}));
        % savefig(fileName)
        saveas(gcf, strcat(fileName, '.png'))
        disp('file saved');
        
        
end



end
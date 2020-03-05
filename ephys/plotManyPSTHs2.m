function plotManyPSTHs2(session, opts)


% settings
s.folder = fullfile(getenv('OBSDATADIR'), 'figures', 'ephys', 'PSTHs');  % folder in which the PSTHs will be saved
s.rows = 4;
s.cols = 5;
s.timeWindow = [-0.05, 0.05];
s.stepPercentiles = [30 70]; % only include steps with durations in between these percentile limits
s.pawNames = {'LH', 'LF', 'RF', 'RH'};
s.pawColors = hsv(4);
s.mainColor = [.8 .4 1];

% reassign settings contained in opts
if exist('opts', 'var'); for i = 1:2:length(opts); s.(opts{i}) = opts{i+1}; end; end


% initializations
sessionFolder = fullfile(getenv('OBSDATADIR'), 'sessions', session);
load(fullfile(sessionFolder, 'runAnalyzed.mat'), ...
        'obsOnTimes', 'obsOffTimes',  'wiskContactFrames', 'frameTimeStamps', ...
        'frameTimeStampsWisk', 'rewardTimes', 'isLightOn', 'touches', 'touchesPerPaw', 'touchClassNames');
disp('finish loading runAnalyzed.mat');        


% load kinData if it exists // otherwise compute kinData    
if exist(fullfile(sessionFolder, 'kinData.mat'))
    load(fullfile(sessionFolder, 'kinData.mat'), 'kinData', 'stanceBins')
else
    [kinData, stanceBins] = getKinematicData(session);
end
disp('finish loading kin data');



% load neural data
load(fullfile(sessionFolder, 'neuralData.mat'));
cellData = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'cellData.csv'));
goodBins = logical([cellData.include]);
cellData = cellData(goodBins,:);
disp('finish loading neural data');


for cellNum = 1:length(unit_ids)
    
    % !!! load cell firing rate and times
    fprintf('%s: plotting cell %i/%i\n', session, cellNum, length(unit_ids))
    figure('name', sprintf('%s - unit %i - 2', session, unit_ids(cellNum)), ...
        'color', 'white', 'MenuBar', 'none', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
    plotInd = 0;
    
    minTime = polyval(openEphysToSpikeMapping, cellData.timeStart(cellNum)*60);
    maxTime = polyval(openEphysToSpikeMapping, cellData.timeEnd(cellNum)*60);
    
    
    % PSTHs center around swing start and stance start timepoints,
    % for each paw
    for paw = 1:4
        % plot PSTH center around swing start point
        plotInd = plotInd + 1; subplot(s.rows, s.cols, plotInd);
        
        swingStartInds = find(diff(~stanceBins(:,paw))==1);
        swingStartTimes = frameTimeStamps(swingStartInds);
        validBins = find(swingStartTimes > minTime & swingStartTimes < maxTime);
        validSwingStartTimes = swingStartTimes(validBins);
        
        swingEndInds = find(diff(stanceBins(:,paw))== 1);
        swingEndTimes = frameTimeStamps(swingEndInds);
        validBins = find(swingEndTimes > minTime & swingEndTimes < maxTime);
        validSwingEndTimes = swingEndTimes(validBins);
        
        % quality check
        if length(validSwingStartTimes) > length(validSwingEndTimes)
            if validSwingStartTimes(end) > validSwingEndTimes(end)
                validSwingStartTimes = validSwingStartTimes(1:end-1);
            else
                validSwingStartTimes = validSwingStartTimes(2:end);
            end
            
        elseif length(validSwingStartTimes) < length(validSwingEndTimes)
            if validSwingEndTimes(1) < validSwingStartTimes(1)
                validSwingEndTimes = validSwingEndTimes(2:end);
            else
                validSwingEndTimes = validSwingEndTimes(1:end-1);
            end
        else
            if validSwingEndTimes(1) < validSwingStartTimes(1)
                validSwingEndTimes = validSwingEndTimes(2:end);
                validSwingStartTimes = validSwingStartTimes(1:end-1);
            end                       
        end
        
        substracted = validSwingEndTimes - validSwingStartTimes;
        if any(substracted < 0)
            fprintf('WARNING: swing start and end timepoints mismatch!');
        end
        
        % discard steps which swing phases are unreasonably long
        if any(substracted > 0.5)
            validSwingEndTimes(find(substracted > 0.5)) = [];
            validSwingStartTimes(find(substracted > 0.5)) = [];
        end
        
        % plot PSTH center around swing start point
        plotPSTH2(session, cellNum, validSwingStartTimes, {'xLims', s.timeWindow, 'colors', s.pawColors(paw,:)});
        xlabel(sprintf('%s: swing start', s.pawNames{paw}))
        
        % plot PSTH center around stance start point
        plotInd = plotInd + 1; subplot(s.rows, s.cols, plotInd);
        plotPSTH2(session, cellNum, validSwingEndTimes, {'xLims', s.timeWindow, 'colors', s.pawColors(paw,:)});
        xlabel(sprintf('%s: stance start', s.pawNames{paw}))
    end
    
     % save
    fileName = fullfile(s.folder, [session 'unit' num2str(unit_ids(cellNum)) '_2']);
%     savefig(fileName)
    saveas(gcf, [fileName '.png'])
    disp('file saved');
    
    
end

    
end
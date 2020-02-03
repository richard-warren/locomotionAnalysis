function plotRawTraces(session, unit_id, varargin)

% settings
s.shadedArea = false; % in case I want to plot shaded area other than obs events.
s.shadedAreaTimeInds = []; % n by 2 matrix, each row is one trial or event, the first colomn is start time ind, the second is stop time ind.
s.yLim = 200; % for plotting
s.xLim = 25; % for plotting
s.neuralActivityColor = [1 0.39 0.43];
s.velocityColor = [0.33 0.73 1];
s.lightOnColor = [1 0.87 0.57];
s.lightOffColor = [0.64 0.24 0.08];
s.rewardLineColor = [0.77 0.42 1];
s.wiskLineColor = [0.15 0.15 0.15];
s.shadedAreaColor = [0.84 0.68 0.94]; % Color for other shaded areas, in case I want to plot shaded area other than obs events.
s.saveFolder = fullfile('Z:\obstacleData\figures\ephys\IndividualTrialOverlayFigs', session);

if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs

% INITIALIZATION
% load runAnalyzed.mat
display(['processing: ', session, ' unit_id' num2str(unit_id)]);
display(['Loading session data: ' session]);
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'frameTimeStamps', 'obsOnTimes', 'obsOffTimes', 'rewardTimes', 'wiskContactFrames', 'isLightOn', ...
    'obsPixPositions', 'frameTimeStampsWisk', 'wiskContactTimes', 'touches', 'touchesPerPaw', 'touchClassNames', ...
    'wheelPositions', 'wheelTimes')
display('finish loading runAnalyzed.mat');

% load kinData.mat

sessionFolder = fullfile(getenv('OBSDATADIR'), 'sessions', session);
if exist(fullfile(sessionFolder, 'kinData.mat'))
    load(fullfile(sessionFolder, 'kinData.mat'), 'kinData', 'stanceBins')
else
    [kinData, stanceBins] = getKinematicData(session);
end
display('finish loading kinData.mat');

% load cellData.csv
cellData = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'cellData.csv'));
goodBins = logical([cellData.include]);
cellData = cellData(goodBins,:);
display('finish loading cellData.csv');

% load the neuralData.mat
% ephysInfo = getSessionEphysInfo(session);
% getVoltage = @(data, channel, inds) double(data.Data.Data(channel,inds))*ephysInfo.bitVolts;
% data = memmapfile(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysInfo.ephysFolder, [ephysInfo.fileNameBase '_CHs.dat']), ...
%     'Format', {'int16', [ephysInfo.channelNum, ephysInfo.smps], 'Data'}, 'Writable', false);
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'neuralData.mat')) ;
neuralTimeStamps = timeStamps; % reference to spike clock
% spkRates = temp.spkRates;
% unit_ids = temp.unit_ids;
fs = round(1/median(diff(neuralTimeStamps)));
display('finish loading neuralData.mat');


% calculate the velodity
vel = getVelocity(wheelPositions, .02, 1000); % should have the same length as wheelTimes
velTimeStamps = wheelTimes;
clear wheelTimes
clear wheelPositions
display('finish calculating velocity');


% restrict the reward trials
display('preprocessing the spkRate data...');

unit_idInd = find(unit_ids == unit_id);
unit_startTime = polyval(openEphysToSpikeMapping, cellData.timeStart(unit_idInd)*60);
unit_endTime = polyval(openEphysToSpikeMapping, cellData.timeEnd(unit_idInd)*60);

rewardTimesSelected = rewardTimes(rewardTimes > unit_startTime & rewardTimes < unit_endTime);
rewardNum = length(rewardTimesSelected);
rewardTimeBinsSelected = [rewardTimesSelected(1:rewardNum-1)-1, rewardTimesSelected(2:rewardNum)+2];


% PLOTTING!!!
display('plotting...');

for i = 1:length(rewardTimeBinsSelected)
    if mod(i, 8) == 1
        figure('name', sprintf('%s_unit%i_plot%i', session, unit_id, i), ...
            'color', 'white', 'MenuBar', 'none', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
        plotInd = 1;
    end
    
    % plot neural activity, from last reward delivery to next reward delivery
    subplot(8, 1, plotInd);
    box off
    trialStartTime = rewardTimeBinsSelected(i, 1);
    trialEndTime = rewardTimeBinsSelected(i, 2);
    neuralActivity = spkRates(unit_idInd, find(neuralTimeStamps > trialStartTime & neuralTimeStamps < trialEndTime));
    velocity = vel(1, find(velTimeStamps > trialStartTime & velTimeStamps < trialEndTime));
    timeDuration = rewardTimeBinsSelected(i, 2) - rewardTimeBinsSelected(i, 1);
    
    x1 = linspace(-1, -1+timeDuration, length(neuralActivity));
    x2 = linspace(-1, -1+timeDuration, length(velocity));
    
   
    yyaxis left
    plot(x1, neuralActivity, '-', 'Color', s.neuralActivityColor, 'LineWidth', 1);
    ylim([0, max(spkRates(unit_idInd, :))+10]);
    xlim([-1, timeDuration-1]);
    ax = gca;
    ax.YColor = s.neuralActivityColor;
    
    % add lines for reward delivery
    line([0, 0], [0, max(spkRates(unit_idInd, :))+10], 'Color', s.rewardLineColor, 'LineWidth', 3);
    line([timeDuration-3, timeDuration-3], [0, max(spkRates(unit_idInd, :))+10], 'Color', s.rewardLineColor, 'LineWidth', 3);
    
    % add lines for wisk contact, add shaded areas for obs events
    trialObsOnTimes = obsOnTimes(find(obsOnTimes > trialStartTime & obsOnTimes < trialEndTime), 1);
    trialObsOffTimes = obsOffTimes(find(obsOffTimes > trialStartTime & obsOffTimes < trialEndTime), 1);
    trialInds = find(obsOnTimes > trialStartTime & obsOnTimes < trialEndTime);
    trialIsLightOn = isLightOn(trialInds, 1);
    trialWiskContactTimes = wiskContactTimes(find(wiskContactTimes > trialStartTime & wiskContactTimes < trialEndTime));
    
    trialWiskContactTimes = trialWiskContactTimes - trialStartTime;
    trialObsOnTimes = trialObsOnTimes - trialStartTime;
    trialObsOffTimes = trialObsOffTimes - trialStartTime;
    
    for j = 1:length(trialObsOnTimes)       
        x = linspace(trialObsOnTimes(j), trialObsOffTimes(j), 100);
        y = repmat(max(spkRates(unit_idInd, :))+10, 1, length(x));
        if trialIsLightOn(j) == 1
            hold on
            h = area(x, y, 'FaceColor', s.lightOnColor);
            set(h, 'LineStyle', 'none');
            h.FaceAlpha = 0.5;
            % set(h, 'Layer', 'bottom');
        else
            hold on
            h = area(x, y, 'FaceColor', s.lightOffColor);
            set(h, 'LineStyle', 'none');
            h.FaceAlpha = 0.5;
            % set(h, 'Layer', 'bottom');
        end
        
        hold on
        line([trialWiskContactTimes(j), trialWiskContactTimes(j)], [0,  max(spkRates(unit_idInd, :))+10], 'Color', s.wiskLineColor, 'LineWidth', 2);
        
    end
    
    
    hold on
    yyaxis right
    plot(x2, velocity, '-', 'Color', s.velocityColor, 'LineWidth', 1);
    box off
    ylim([0 max(velocity)+0.1]);
    ax = gca;
    ax.YColor = s.velocityColor;
    % title(['trail ' num2str(i)]);
    xlabel(['time ' '- trial ' num2str(find(rewardTimes == trialStartTime+1))]); hold on;          
   
    plotInd = plotInd + 1;
    
    if plotInd == 9
        savefig(fullfile(s.saveFolder, [session '_unit' num2str(unit_id) '_plot0' num2str(floor(i/8))]));
        saveas(gcf, fullfile(s.saveFolder, [session '_unit' num2str(unit_id) '_plot0' num2str(floor(i/8)) '.png']));
    end
    
    
end 

end    
  
    
    
    
    







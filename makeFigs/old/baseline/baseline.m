function baseline(mouse)

% BASELINE ANALYSIS

% iterates over all sessions of given mouse with experiment name 'baseline'
% plots mean speed as a function of position with respect to reward
% also plots bar graph of average median speed per session to determine whether they are running well enough to include in subsequent experiments
% trial median speeds are computed in space between vertical lines on the plot (positRange)
% only middle trials are included, defined by trialRange



% user settings
dataDir = [getenv('OBSDATADIR') 'sessions\'];
resultsDir = [getenv('OBSDATADIR') 'figures\'];
trialRange = [.05 1]; % only include trials in the middle between these two limits
rewardRotations = 9.1;
positRange = [1 7]; % units: wheel rotations // only compute trial median velocity within these wheel positions on a per-trial basis 
wheelDiam = 0.1905; % m
ylims = [0 .7];



% initializations
sessionInfo = readtable([dataDir 'sessionInfo.xlsx']);

maxPosit = pi * wheelDiam * rewardRotations;
positRangeMeters = pi * wheelDiam * positRange;
sessionInfo.include(isnan(sessionInfo.include)) = 0;
sessionInds = strcmp(sessionInfo.mouse, mouse) &...
              strcmp(sessionInfo.experiment, 'baseline') &...
              sessionInfo.include;
sessions = sessionInfo.session(sessionInds);

cmap = winter(length(sessions));
figure('name', mouse);
subplot(1,2,2); bar(nan(1,length(sessions))); hold on % ghost bar plot to get our axis labels


% plot sessions means
for i = 1:length(sessions)
    
    % load session data
    load([dataDir sessions{i} '\runAnalyzed.mat'],...
        'wheelPositions', 'wheelTimes', 'rewardTimes', 'targetFs')
    
    positsInterp = 0 : (1/targetFs) : maxPosit;
    
    
    % limit to middle trials only
    trialLims = round(trialRange * length(rewardTimes));
    trialLims = max(trialLims, 1); % ensure no 0 indices
    rewardTimes = rewardTimes(trialLims(1):trialLims(2));
    

    % compute velocity
    vel = getVelocity(wheelPositions, .5, targetFs);

    
    % get per trial velocity and positions (cell arrays with one trial per entry)
    vel = splitByRewards(vel, wheelTimes, rewardTimes, false);
    posits = splitByRewards(wheelPositions, wheelTimes, rewardTimes, true);

    
    % interpolate velocities over evenly spaced positional values
    velInterp = nan(length(rewardTimes), length(positsInterp));

    for j = 1:length(posits)

        % remove duplicate positional values
        [posits{j}, uniqueInds] = unique(posits{j}, 'stable');
        vel{j} = vel{j}(uniqueInds);

        % interpolate
        velInterp(j,:) = interp1(posits{j}, vel{j}, positsInterp, 'linear');

    end
    
    
    % compute and plot session average
    sessionMean = nanmean(velInterp, 1);
    subplot(1,2,1)
    plot(positsInterp, sessionMean, 'color', cmap(i,:), 'linewidth', 3)
    hold on
    
    % compute and plot median velocity
    subplot(1,2,2)
    middlePositInds = (positsInterp > positRangeMeters(1)) & (positsInterp < positRangeMeters(2));
    trialMedians = nanmedian(velInterp(:, middlePositInds), 2);
    sessionMed = nanmedian(trialMedians);
    bar(i, sessionMed, 'facecolor', cmap(i,:)); hold on
end



% pimp figure
pimpFig; set(gcf, 'menubar', 'none', 'position', [.2 .4 .6 .4])

subplot(1,2,1); set(gca, 'ylim', ylims, 'xlim', [0 maxPosit])
for i=1:2
    line([positRangeMeters(i) positRangeMeters(i)], ylims, 'color', [.2 .2 .2]);
end
xlabel('distance travelled (m)', 'fontweight', 'bold')
ylabel('velocity (m/s)', 'fontweight', 'bold')

subplot(1,2,2);
set(gca, 'ylim', ylims, 'xlim', [.25 length(sessions)+.75])
xlabel('session #', 'fontweight', 'bold')
ylabel('velocity (m/s)', 'fontweight', 'bold')



% save figure
savefig([resultsDir mouse 'Baseline.fig'])



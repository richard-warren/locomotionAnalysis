
% PLOT BASELINE DATA



% user settings
dataDir = 'C:\Users\Rick\Google Drive\columbia\obstacleData\sessions\';
mouse = 'run4';
trialRange = [.2 .8]; % only include trials in the middle between these two limits
rewardRotations = 8;
wheelDiam = 0.1905; % m
ylims = [0 .5];



% initializations
sessionInfo = readtable([dataDir 'sessionInfo.xlsx']);

maxPosit = pi * wheelDiam * rewardRotations;
sessionInds = strcmp(sessionInfo.mouse, mouse) &...
              strcmp(sessionInfo.experiment, 'baseline') &...
              sessionInfo.include;
sessions = sessionInfo.session(sessionInds);

cmap = copper(length(sessions));
figure;



% plot sessions means
for i = 1:length(sessions)
    
    % load session data
    load([dataDir sessions{i} '\runAnalyzed.mat'],...
        'wheelPositions', 'wheelTimes', 'rewardTimes', 'targetFs')
    
    positsInterp = 0 : (1/targetFs) : maxPosit;
    
    
    % trim first and last rewards
    lims = round(trialRange * length(rewardTimes));
    rewardTimes = rewardTimes(lims(1):lims(2));
    

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
    plot(positsInterp, sessionMean, 'color', cmap(i,:), 'linewidth', 2)
    hold on
end

pimpFig;
set(gca, 'ylim', ylims)


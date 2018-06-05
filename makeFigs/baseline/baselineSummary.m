function baselineSummary(mice)

% BASELINE ANALYSIS

% iterates over all sessions for all mice with experiment name 'baseline'
% plots mean speed as a function of position with respect to reward for each mouse, and average across mice
% also plots 1D scatter graph of mean speed per mouse (across sessions)
% mean speeds are computed in space between vertical lines on the plot (positRange)
% only middle trials are included, defined by trialRange



% user settings
dataDir = [getenv('OBSDATADIR') 'sessions\'];
resultsDir = [getenv('OBSDATADIR') 'figures\'];
trialRange = [.05 .95]; % only include trials in the middle between these two limits
rewardRotations = 9.1;
positRange = [2 7]; % units: wheel rotations // only compute trial median velocity within these wheel positions on a per-trial basis 
wheelDiam = 0.1905; % m
yLims = [0 .8];
smoothing = .5; % window for mean smoothing (s)
yTicks = 0:.2:yLims(2);



% initializations
sessionInfo = readtable([dataDir 'sessionInfo.xlsx']);

maxPosit = pi * wheelDiam * rewardRotations;
positRangeMeters = pi * wheelDiam * positRange;
% sessionInfo.include(isnan(sessionInfo.include)) = 0;
sessionBins = ismember(sessionInfo.mouse, mice) &...
              strcmp(sessionInfo.experiment, 'baseline') &...
              sessionInfo.include;
sessions = sessionInfo(sessionBins, :);

data = struct(); % initialize storage variable
cmap = hsv(length(mice));


% collect data
for i = 1:size(sessions,1)
    
    % load session data
    load([dataDir sessions.session{i} '\runAnalyzed.mat'],...
        'wheelPositions', 'wheelTimes', 'rewardTimes', 'targetFs')
    
    positsInterp = 0 : (1/targetFs) : maxPosit;
    smoothSmps = targetFs * smoothing;
    
    
    % limit to middle trials only
    trialLims = round(trialRange * length(rewardTimes));
    trialLims = max(trialLims, 1); % ensure no 0 indices
    rewardTimes = rewardTimes(trialLims(1):trialLims(2));
    

    % compute velocity
    velContinuous = getVelocity(wheelPositions, .5, targetFs);

    
    % get per trial velocity and positions (cell arrays with one trial per entry)
    vel = splitByRewards(velContinuous, wheelTimes, rewardTimes, false);
    pos = splitByRewards(wheelPositions, wheelTimes, rewardTimes, true);

    
    % interpolate velocities over evenly spaced positional values
    velInterp = nan(length(rewardTimes), length(positsInterp));

    for j = 1:length(pos)

        % remove duplicate positional values
        [pos{j}, uniqueInds] = unique(pos{j}, 'stable');
        vel{j} = vel{j}(uniqueInds);

        % interpolate
        velInterp(j,:) = interp1(pos{j}, vel{j}, positsInterp, 'linear');

    end
    
    
    % store values
    data(i).mouse = sessions.mouse{i};
    data(i).vel = nanmean(velInterp, 1);
    data(i).vel = smooth(data(i).vel, smoothSmps);
    middlePositInds = (positsInterp > positRangeMeters(1)) & (positsInterp < positRangeMeters(2));
    data(i).meanVel = nanmean( nanmean(velInterp(:, middlePositInds), 2) ); % !!! why is nanmean necessary here!!!
    
end


% plot mouse means


% compute and plot mouse means

figure('name', 'baselineSummary', 'menubar', 'none', 'units', 'pixels', 'position', [500 500 900 400], 'color', [1 1 1]);

allVels = reshape([data.vel], length([data(1).vel]), length(data))';
mouseMeanTraces = nan(length(mice), size(allVels,2));
mouseMedVels = nan(length(mice), 1);
jitterPos = linspace(-.25,.25,length(mice));
jitterPos = jitterPos(randperm(length(jitterPos)));

for i = 1:length(mice)
    
    % get mouse means across sessions
    bins = strcmp(mice{i}, {data.mouse});
    mouseMeanTraces(i,:) = nanmean(allVels(bins,:),1);
    mouseMedVels(i) = mean([data(bins).meanVel]);
    
    % trace
    subplot(1,3,1:2)
    plot(positsInterp, mouseMeanTraces(i,:), 'linewidth', 1.5, 'color', cmap(i,:)); hold on
    
    % scatter
    subplot(1,3,3)
    scatter(jitterPos(i), mouseMedVels(i), 75, cmap(i,:), 'filled'); hold on
    
end


% plot means across mice

subplot(1,3,1:2)
plot(positsInterp, nanmean(mouseMeanTraces,1), 'color', [0 0 0], 'linewidth', 4)

subplot(1,3,3)
meanVel = mean(mouseMedVels);
line([-.5 .5], [meanVel meanVel], 'color', [0 0 0], 'linewidth', 3)

legend(mice, 'location', 'southeast')
legend('boxoff')


% pimp fig

% velocity traces
subplot(1,3,1:2); set(gca, 'ylim', yLims, 'xlim', [0 maxPosit], 'box', 'off', 'xtick', 0:1:5, 'ytick', yTicks)

% add lines indicating position range for mean calculation
for i = 1:2
    line([positRangeMeters(i) positRangeMeters(i)], [0 .05], 'linewidth', 2, 'color', get(gca, 'xcolor'));
end

xlabel('distance travelled (m)', 'fontweight', 'bold')
ylabel('velocity (m/s)', 'fontweight', 'bold')

% scatter plot
subplot(1,3,3); set(gca, 'ylim', yLims, 'xlim', [-1.5 3], 'xcolor', 'none', 'ytick', yTicks)
ylabel('velocity (m/s)', 'fontweight', 'bold')


% save figure
savefig([resultsDir 'baselineSummary.fig'])



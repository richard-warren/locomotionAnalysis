% OBSTACLE AVOIDANCE ANALYSIS

% note: if i ever add obstacle on-position jitter this code won't work // i will have to incorporate the wheel rotary encoder, or track backwords from the position at which the obs turn OFF, which remains constant

% user settings
dataDir = 'C:\Users\Rick\Google Drive\columbia\obstacleData\sessions\';
mouse = 'run3';

obsPrePost = [.1 .5]; % meters // plot this much before and after the obstacle turns on
obsOnPos = 0.4236; % meters // distance for which obstacle is lit up // this was hack-ily determined emperically in matlab... knowing this value ahead of time makes the coding much simpler though...
posRes = .001; % resolution of x axis, in meters
ylims = [.1 .6]; % m/s
trialRange = [0 .8]; % only include trials #s between these limits, so performance metrices are not dragged down when the animal is warming up or sated near the end



% initializations
sessionInfo = readtable([dataDir 'sessionInfo.xlsx']);

sessionInds = strcmp(sessionInfo.mouse, mouse) &...
              strcmp(sessionInfo.experiment, 'obsTest') &...
              sessionInfo.include;
sessions = sessionInfo.session(sessionInds);

posInterp = -obsPrePost(1) : posRes : (obsOnPos+obsPrePost(2)); % velocities will be interpolated across this grid of positional values

cmap = winter(length(sessions));
figure;
subplot(2,2,2); bar(nan(1,length(sessions))); hold on % ghost bar plot to get our axis labels



% iterate over sessions
for i = 1:length(sessions)

    % load session data
    load([dataDir sessions{i} '\run.mat'], 'ObsLight', 'touch');
    load([dataDir sessions{i} '\runAnalyzed.mat'],...
            'wheelPositions', 'wheelTimes', 'rewardTimes', 'targetFs');
    obsOnTimes = ObsLight.times(logical(ObsLight.level)); % important: assumes first event is HIGH... not sure how this will behave otherwise...
    
    % limit to middle trials only
    trialLims = round(trialRange * length(obsOnTimes));
    trialLims = max(trialLims, 1); % ensure no 0 indices
    obsOnTimes = obsOnTimes(trialLims(1):trialLims(2));
    
    % get touch times // !!! i should really debounce these signals so the duration of touch is interpretable... currently duration is not interpretable
    touchTimes = touch.times(logical([0; diff(touch.values>3)==1]));

    % compute velocity
    vel = getVelocity(wheelPositions, .5, targetFs);

    
    
    % iterate over all trials, computing velocity vs. position
    trialVels = cell(length(obsOnTimes));
    sessionVels = nan(length(obsOnTimes), length(posInterp));
    obsAvoided = nan(1,length(obsOnTimes));

    for j = 1:length(obsOnTimes)

        % locate trial
        onInd  = find(wheelTimes > obsOnTimes(j), 1, 'first');    
        onPos = wheelPositions(onInd);

        % get trial positions and velocities
        trialInds = (wheelPositions > onPos-obsPrePost(1)) & (wheelPositions < onPos+obsOnPos+obsPrePost(2));
        trialPos = wheelPositions(trialInds);
        trialPos = trialPos - trialPos(1); % baseline normalize
        trialVel = vel(trialInds);

        % remove duplicate positional values
        [trialPos, uniqueInds] = unique(trialPos, 'stable');
        trialVel = trialVel(uniqueInds);

        % interpolate velocities across positional grid
        trialVelInterp = interp1(trialPos, trialVel, posInterp, 'linear', 'extrap');

        % store results
        sessionVels(j,:) = trialVelInterp;
        
        % determine whether obstacle was avoided
        onTime =  wheelTimes(find(wheelPositions>onPos,1,'first'));
        offTime = wheelTimes(find(wheelPositions<(onPos+obsOnPos),1,'last'));
        obsAvoided(j) = ~any(touchTimes>onTime & touchTimes<offTime);
    end
    
    % plot session mean velocity
    subplot(2,2,1)
    plot(posInterp, mean(sessionVels,1), 'color', cmap(i,:), 'linewidth', 2); hold on
    
    % plot session success rate
    subplot(2,2,2)
    bar(i, mean(obsAvoided), 'facecolor', cmap(i,:)); hold on
    
    % plot successful trial velocity means
    subplot(2,2,3)
    plot(posInterp, mean(sessionVels(logical(obsAvoided),:),1), 'color', cmap(i,:), 'linewidth', 2); hold on
    
    % plot unsuccessful trial velocity means
    subplot(2,2,4)
    plot(posInterp, mean(sessionVels(~logical(obsAvoided),:),1), 'color', cmap(i,:), 'linewidth', 2); hold on
    
end



% pimp out figure
pimpFig; set(gcf, 'menubar', 'none', 'position', [.2 .1 .6 .8])

% pimp out all velocity figures
titles = {'all trials', 'successful trials', 'unsuccessful trials'}; % the empty slot is sort of a hack
velPlots = [1,3,4];

for i = 1:length(velPlots)
    
    subplot(2,2,velPlots(i))
    set(gca, 'xlim', [-obsPrePost(1) (obsOnPos+obsPrePost(2))], 'ylim', ylims)
    
    title(titles(i));
    xlabel('position (m)', 'fontweight', 'bold')
    ylabel('velocity (m/s)', 'fontweight', 'bold')
    
    line([0 0], ylims, 'color', [0 0 0])
    line([obsOnPos obsOnPos], ylims, 'color', [0 0 0])
    
end

% add titles
subplot

% pimp out bar graph
subplot(2,2,2)
set(gca, 'ylim', [0 1], 'xlim', [.25 length(sessions)+.75])
title('success rate')
xlabel('session #', 'fontweight', 'bold')
ylabel('fraction avoided', 'fontweight', 'bold')



% save fig
savefig(['obsAvoidance\figs\' mouse '.fig'])




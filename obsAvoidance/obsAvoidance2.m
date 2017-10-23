function obsAvoidance2(mouse, expName)

% OBSTACLE AVOIDANCE ANALYSIS

% user settings
% dataDir = 'C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\sessions\';
dataDir = 'C:\Users\Rick\Google Drive\columbia\obstacleData\sessions\';

obsPrePost = [.6 .25]; % plot this much before and after the obstacle reaches the mouse
posRes = .001; % resolution of x axis, in meters
touchPosRes = .0001;
ylims = [.1 .6]; % m/s
trialRange = [0 .8]; % only include trials #s between these limits, so performance metrices are not dragged down when the animal is warming up or sated near the end
obsPos = .382; % m, position at which obstacle is in the middle of the frame // use getFrameTimes function to determine this value
frameEdges = [.336 .444]; % (m)
sig = .005; % sigma for gaussian kernal


% initializations
sessionInfo = readtable([dataDir 'sessionInfo.xlsx']);

sessionInds = strcmp(sessionInfo.mouse, mouse) &...
              strcmp(sessionInfo.experiment, expName) &...
              sessionInfo.include;
sessions = sessionInfo.session(sessionInds);

posInterp = -obsPrePost(1) : posRes : obsPrePost(2); % velocities will be interpolated across this grid of positional values
touchPosInterp = frameEdges(1): touchPosRes : frameEdges(2);

gausKernel = arrayfun(@(x) (1/(sig*sqrt(2*pi))) * exp(-.5*(x/sig)^2), -sig*5:1/(1/touchPosRes):sig*5);
gausKernel = gausKernel/sum(gausKernel);

cmap = winter(length(sessions));
figure;
% subplot(2,2,2); bar(nan(1,length(sessions))); hold on % ghost bar plot to get our axis labels



% iterate over sessions
for i = 1:length(sessions)

    % load session data
    load([dataDir sessions{i} '\run.mat'], 'touch');
    load([dataDir sessions{i} '\runAnalyzed.mat'],...
            'wheelPositions', 'wheelTimes',...
            'obsPositions', 'obsTimes',...
            'obsOnTimes', 'obsOffTimes',...
            'rewardTimes', 'targetFs');
    
    obsPositions = fixObsPositions(obsPositions, obsTimes, obsOnTimes);
    
    % limit to middle trials only
    trialLims = round(trialRange * length(obsOnTimes));
    trialLims = max(trialLims, 1); % ensure no 0 indices
    obsOnTimes = obsOnTimes(trialLims(1):trialLims(2));
    obsOffTimes = obsOffTimes(trialLims(1):trialLims(2));
    
    % get touch times
    touchTimes = touch.times(logical([0; diff(touch.values>3)==1]));

    % compute velocity
    vel = getVelocity(wheelPositions, .5, targetFs);

    
    
    % iterate over all trials
    sessionVels = nan(length(obsOnTimes), length(posInterp));
    obsAvoided = nan(1,length(obsOnTimes));
    touchPositions = [];
    firstTouchTimes = [];
    
    obsOnPositions = nan(1,length(obsOnTimes)); % record positions at which obstacle turns on (this may be jittered)

    for j = 1:length(obsOnTimes)

        % VELOCITY VS POSITION
        
        % locate trial
        obsOnPos = obsPositions( find(obsTimes >= obsOnTimes(j), 1, 'first') );
        obsTime  = obsTimes(find( obsTimes >= obsOnTimes(j) & obsPositions >= obsPos, 1, 'first')); % time at which obstacle reaches obsPos
        obsWheelPos = wheelPositions(find(wheelTimes>=obsTime, 1, 'first')); % position of wheel at moment obstacle reaches obsPos

        % get trial positions and velocities
        trialInds = (wheelPositions > obsWheelPos-obsPrePost(1)) & (wheelPositions < obsWheelPos+obsPrePost(2));
        trialPos = wheelPositions(trialInds);
        trialPos = trialPos - obsWheelPos; % normalize s.t. 0 corresponds to the position at which the obstacle is directly over the wheel
        trialVel = vel(trialInds);

        % remove duplicate positional values
        [trialPos, uniqueInds] = unique(trialPos, 'stable');
        trialVel = trialVel(uniqueInds);

        % interpolate velocities across positional grid
        trialVelInterp = interp1(trialPos, trialVel, posInterp, 'linear');

        % store results
        sessionVels(j,:) = trialVelInterp;
        
        % find whether and where obstacle was touched
        touchInd = find(touchTimes>obsOnTimes(j) & touchTimes<obsOffTimes(j), 1, 'first');
        obsAvoided(j) = ~any(touchInd);
        touchPos = obsPositions(find(obsTimes>=touchTimes(touchInd), 1, 'first'));
        if ~isempty(touchPos)
            touchPositions(end+1) = touchPos;
            firstTouchTimes(end+1) = touchTimes(touchInd);
        end
        
        % record position at which obstacle turned on
        obsOnPositions(j) = obsPos - obsOnPos;
    end
    
    
    % compute touch probability
    touchCounts = histcounts(touchPositions, touchPosInterp);
    touchCounts = touchCounts / length(obsOnTimes);
    touchProb = conv(touchCounts, gausKernel, 'same');
    
    % compute avg frame at moment of all first touches
    avgTouchFrame = getAvgFrameAtTimes(sessions{i}, 'Bot', firstTouchTimes);
    
    % plot touch probability
    subplot(2,1,1)
    touchPosCenters = touchPosInterp(1:end-1) - .5*touchPosRes;
    plot(touchPosCenters, touchProb, 'color', cmap(i,:), 'linewidth', 2); hold on;
    

    
%     obsOnMean = mean(obsOnPositions);
    
%     % plot session mean velocity
%     subplot(2,2,1)
%     plot(posInterp, nanmean(sessionVels,1), 'color', cmap(i,:), 'linewidth', 2); hold on
%     line(-[obsOnMean obsOnMean], ylims, 'color', cmap(1,:));
%     
%     % plot session success rate
%     subplot(2,2,2)
%     bar(i, mean(obsAvoided), 'facecolor', cmap(i,:)); hold on
%     
%     % plot successful trial velocity means
%     subplot(2,2,3)
%     plot(posInterp, nanmean(sessionVels(logical(obsAvoided),:),1), 'color', cmap(i,:), 'linewidth', 2); hold on
%     line(-[obsOnMean obsOnMean], ylims, 'color', cmap(1,:));
%     
%     % plot unsuccessful trial velocity means
%     subplot(2,2,4)
%     plot(posInterp, nanmean(sessionVels(~logical(obsAvoided),:),1), 'color', cmap(i,:), 'linewidth', 2); hold on
%     line(-[obsOnMean obsOnMean], ylims, 'color', cmap(1,:));
end

% pimp out touch probability
subplot(2,1,1)
set(gca, 'xdir', 'reverse', 'xlim', [touchPosCenters(1) touchPosCenters(end)]);

% show average touch frame
subplot(2,1,2)
frame = getAvgFrameAtObsLocation('171020_002', obsPos);
scaling = diff(frameEdges) / size(frame,2);
imagesc(frameEdges, [0 size(frame,1)*scaling], fliplr(imadjust(avgTouchFrame)));
set(gca, 'xdir', 'reverse', 'xlim', [touchPosCenters(1) touchPosCenters(end)],...
    'xtick', [], 'ytick', [])

% pimp out figure
pimpFig; set(gcf, 'menubar', 'none', 'position', [.2 .1 .6 .8])
keyboard
% 
% % pimp out all velocity figures
% titles = {'all trials', 'successful trials', 'unsuccessful trials'}; % the empty slot is sort of a hack
% velPlots = [1,3,4];
% 
% for i = 1:length(velPlots)
%     
%     subplot(2,2,velPlots(i))
%     set(gca, 'xlim', [-obsPrePost(1) obsPrePost(2)], 'ylim', ylims)
%     
%     title(titles(i));
%     xlabel('position (m)', 'fontweight', 'bold')
%     ylabel('velocity (m/s)', 'fontweight', 'bold')
%     
%     line([0 0], ylims, 'color', cmap(1,:))
%     
% end
% 
% % pimp out bar graph
% subplot(2,2,2)
% set(gca, 'ylim', [0 1], 'xlim', [.25 length(sessions)+.75])
% title('success rate')
% xlabel('session #', 'fontweight', 'bold')
% ylabel('fraction avoided', 'fontweight', 'bold')
% 


% save fig
savefig(['obsAvoidance\figs\' mouse '.fig'])

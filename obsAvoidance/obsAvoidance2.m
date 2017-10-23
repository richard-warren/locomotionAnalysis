function obsAvoidance2(mouse, expName)

% OBSTACLE AVOIDANCE ANALYSIS
%
% input         mouse:      name of mouse to analyze
%               expName:    string or cell array of experiments to include in analysis

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
sig = .0025; % sigma for gaussian kernal
probYlims = [0 .008];


% initializations
sessionInfo = readtable([dataDir 'sessionInfo.xlsx']);

sessionInds = strcmp(sessionInfo.mouse, mouse) &...
              cellfun(@(x) any(strcmp(x, expName)), sessionInfo.experiment) &...
              sessionInfo.include;
sessions = sessionInfo.session(sessionInds);

posInterp = -obsPrePost(1) : posRes : obsPrePost(2); % velocities will be interpolated across this grid of positional values
touchPosInterp = frameEdges(1): touchPosRes : frameEdges(2);

gausKernel = arrayfun(@(x) (1/(sig*sqrt(2*pi))) * exp(-.5*(x/sig)^2), -sig*5:1/(1/touchPosRes):sig*5);
gausKernel = gausKernel/sum(gausKernel);

cmap = winter(length(sessions));
figure;
subplot(2,2,4); bar(nan(1,length(sessions))); hold on % ghost bar plot to get our axis labels



% iterate over sessions
for i = 1:length(sessions)

    % load session data
    load([dataDir sessions{i} '\runAnalyzed.mat'],...
            'wheelPositions', 'wheelTimes',...
            'obsPositions', 'obsTimes',...
            'obsOnTimes', 'obsOffTimes',...
            'touchOnTimes',...
            'rewardTimes', 'targetFs');
    
    obsPositions = fixObsPositions(obsPositions, obsTimes, obsOnTimes);
    
    % get touch positions and ensure all touches fall within frame
    touchPositions = interp1(obsTimes, obsPositions, touchOnTimes, 'linear');
    validInds = touchPositions>frameEdges(1) & touchPositions<frameEdges(2);
    touchOnTimes = touchOnTimes(validInds);
    touchPositions = touchPositions(validInds);
    
    
    % limit to middle trials only
    trialLims = round(trialRange * length(obsOnTimes));
    trialLims = max(trialLims, 1); % ensure no 0 indices
    obsOnTimes = obsOnTimes(trialLims(1):trialLims(2));
    obsOffTimes = obsOffTimes(trialLims(1):trialLims(2));
    

    % compute velocity
    vel = getVelocity(wheelPositions, .5, targetFs);
    
    
    % iterate over all trials
    sessionVels = nan(length(obsOnTimes), length(posInterp));
    obsAvoided = nan(1,length(obsOnTimes));
    
    obsOnPositions = nan(1,length(obsOnTimes)); % record positions at which obstacle turns on (this may be jittered)

    for j = 1:length(obsOnTimes)
        
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
        
        % find whether obstacle was avoided
        obsAvoided(j) = ~any(touchOnTimes>obsOnTimes(j) & touchOnTimes<obsOffTimes(j));
        
        % record position at which obstacle turned on
        obsOnPositions(j) = obsPos - obsOnPos;
    end
    
    
    % compute touch probability
    touchCounts = histcounts(touchPositions, touchPosInterp);
    touchCounts = touchCounts / length(obsOnTimes);
    touchProb = conv(touchCounts, gausKernel, 'same');
    
    % compute avg frame at moment of all first touches (only for most recent session)
    if i==length(sessions)
        avgTouchFrame = getAvgFrameAtTimes(sessions{i}, 'Bot', touchOnTimes);
    end
    
    % plot touch probability
    subplot(2,2,1)
    touchPosCenters = touchPosInterp(1:end-1) - .5*touchPosRes;
    plot(touchPosCenters - obsPos, touchProb, 'color', cmap(i,:), 'linewidth', 2); hold on;
    

    % get mean obstacle start position
    obsOnMean = mean(obsOnPositions);
    
    % plot session mean velocity
    subplot(2,2,2)
    plot(posInterp, nanmean(sessionVels,1), 'color', cmap(i,:), 'linewidth', 2); hold on
    line(-[obsOnMean obsOnMean], ylims, 'color', cmap(i,:), 'linewidth', 2);
    
    % plot session success rate
    subplot(2,2,4)
    bar(i, mean(obsAvoided), 'facecolor', cmap(i,:)); hold on

end




% pimp out figure
pimpFig;
set(gcf, 'menubar', 'none',...
         'units', 'inches',...
         'position', [4 1.5 10 6.5])


% pimp out touch probability
subplot(2,2,1)
set(gca, 'xdir', 'reverse', 'xlim', [touchPosCenters(1) touchPosCenters(end)] - obsPos, 'ylim', probYlims);
title('touch probability')
xlabel('\itposition (m)')
ylabel('\ittouch probability')


% plot average touch frame
subplot(2,2,3)
frame = getAvgFrameAtObsLocation('171020_002', obsPos);
scaling = diff(frameEdges) / size(frame,2);
imagesc(frameEdges, [0 size(frame,1)*scaling], fliplr(imadjust(avgTouchFrame)));
set(gca, 'xdir', 'reverse', 'xlim', [touchPosCenters(1) touchPosCenters(end)],...
    'xtick', [], 'ytick', [])

 
% pimp out velocity figure
subplot(2,2,2)
set(gca, 'xlim', [-obsPrePost(1) obsPrePost(2)], 'ylim', ylims)

title('velocity');
xlabel('\itposition (m)')
ylabel('\itvelocity (m/s)')
y1 = frameEdges(1)-obsPos;
y2 = frameEdges(2)-obsPos;
line([y1 y1], ylims, 'color', cmap(1,:), 'linewidth', 2)
line([y2 y2], ylims, 'color', cmap(1,:), 'linewidth', 2)
    

% pimp out bar graph
subplot(2,2,4)
set(gca, 'ylim', [0 1], 'xlim', [.25 length(sessions)+.75])
title('success rate')
xlabel('\itsession #')
ylabel('\itfraction avoided')


% save fig
savefig(['obsAvoidance\figs\' mouse '.fig'])

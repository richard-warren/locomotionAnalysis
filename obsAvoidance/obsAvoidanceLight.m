function obsAvoidanceLight(mouse, expName)


% OBSTACLE AVOIDANCE LIGHT
%
% input         mouse:      name of mouse to analyze
%               expName:    string or cell array of experiments to include in analysis




% user settings
% dataDir = 'C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\sessions\';
dataDir = 'C:\Users\Rick\Google Drive\columbia\obstacleData\sessions\';

obsPrePost = [.6 .25]; % plot this much before and after the obstacle reaches the mouse
posRes = .001; % resolution of x axis, in meters
touchPosRes = .0001;
ylims = [.1 .7]; % m/s
trialRange = [0 .8]; % only include trials #s between these limits, so performance metrices are not dragged down when the animal is warming up or sated near the end
obsPos = .382; % m, position at which obstacle is in the middle of the frame // use getFrameTimes function to determine this value
frameEdges = [.336 .444]; % (m)
sig = .0025; % sigma for gaussian kernal
probYlims = [0 .008];
minTouchTime = 0; % only touches count that are >= minTouchTime


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

cmap{1} = winter(length(sessions));
cmap{2} = copper(length(sessions));
figure('name', mouse);
subplot(2,3,3); bar(nan(1,length(sessions))); hold on % ghost bar plot to get our axis labels
subplot(2,3,6); bar(nan(1,length(sessions))); hold on % ghost bar plot to get our axis labels
labels = {' (light on)', ' (light off)'};



% iterate over sessions
for i = 1:length(sessions)

    % load session data
    load([dataDir sessions{i} '\runAnalyzed.mat'],...
            'wheelPositions', 'wheelTimes',...
            'obsPositions', 'obsTimes',...
            'obsOnTimes', 'obsOffTimes',...
            'obsLightOnTimes', 'obsLightOffTimes',...
            'touchOnTimes', 'touchOffTimes',...
            'rewardTimes', 'targetFs');
    
    obsPositions = fixObsPositions(obsPositions, obsTimes, obsOnTimes);
    
    
    % remove brief touches
    validLengthInds = (touchOffTimes - touchOnTimes) >= minTouchTime;
    touchOnTimes = touchOnTimes(validLengthInds);
    touchOffTimes = touchOffTimes(validLengthInds);
    
    
    % get touch positions and ensure all touches fall within frame
    touchPositions = interp1(obsTimes, obsPositions, touchOnTimes, 'linear');
    validPosInds = touchPositions>frameEdges(1) & touchPositions<frameEdges(2);
    touchOnTimes = touchOnTimes(validPosInds);
    touchOffTimes = touchOffTimes(validPosInds);
    touchPositions = touchPositions(validPosInds);
    
    
    % limit to middle trials only
    trialLims = round(trialRange * length(obsOnTimes));
    trialLims = max(trialLims, 1); % ensure no 0 indices
    obsOnTimes = obsOnTimes(trialLims(1):trialLims(2));
    obsOffTimes = obsOffTimes(trialLims(1):trialLims(2));
    
    trialLims = round(trialRange * length(obsLightOnTimes));
    trialLims = max(trialLims, 1); % ensure no 0 indices
    obsLightOnTimes = obsLightOnTimes(trialLims(1):trialLims(2));
    obsLightOffTimes = obsLightOffTimes(trialLims(1):trialLims(2));
    
    
    % compute velocity
    vel = getVelocity(wheelPositions, .5, targetFs);
    
    
    % iterate over all trials
    sessionVels = nan(length(obsOnTimes), length(posInterp));
    obsAvoided = nan(1,length(obsOnTimes));
    obsLightOn = nan(1,length(obsOnTimes));
    obsTouchLightOn = nan(1,length(touchPositions));
    
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
        
        % find whether light was on
        obsLightOn(j) = min(abs(obsOnTimes(j) - obsLightOnTimes)) < .5;
        
        % record whether light was on for touches in trial
        obsTouchLightOn(touchOnTimes>=obsOnTimes(j) & touchOnTimes<=obsOffTimes(j)) = obsLightOn(j);
        
        % record position at which obstacle turned on
        obsOnPositions(j) = obsPos - obsOnPos;
    end
    
    
    for k = 1:2
        
        % set inds to plot light on (k=1) or light off (k=2) trials
        if k==1
            inds = obsLightOn;
            touchInds = (obsTouchLightOn==1);
        else
            inds = ~obsLightOn;
            touchInds = (obsTouchLightOn==0);
        end
        
        % compute touch probability
        touchCounts = histcounts(touchPositions(touchInds), touchPosInterp);
        touchCounts = touchCounts / sum(inds);
        touchProb = conv(touchCounts, gausKernel, 'same');

        
        % plot touch probability
        subplot(2,3,1 + (k-1)*3)
        touchPosCenters = touchPosInterp(1:end-1) - .5*touchPosRes;
        plot(touchPosCenters - obsPos, touchProb, 'color', cmap{k}(i,:), 'linewidth', 2); hold on;


        % get mean obstacle start position
        obsOnMean = mean(obsOnPositions);

        % plot session mean velocity
        subplot(2,3,2 + (k-1)*3)
        plot(posInterp, nanmean(sessionVels(logical(inds),:), 1), 'color', cmap{k}(i,:), 'linewidth', 2); hold on
        line(-[obsOnMean obsOnMean], ylims, 'color', cmap{k}(i,:), 'linewidth', 2);

        % plot session success rate
        subplot(2,3,3 + (k-1)*3)
        bar(i, mean(obsAvoided(logical(inds))), 'facecolor', cmap{k}(i,:)); hold on
    end
end




% pimp out figure
pimpFig;
set(gcf, 'menubar', 'none',...
         'units', 'inches',...
         'position', [4 1.5 11 6.5]);


for k = 1:2
    
    % pimp out touch probability
    subplot(2,3,1 + (k-1)*3)
    set(gca, 'xdir', 'reverse', 'xlim', [touchPosCenters(1) touchPosCenters(end)] - obsPos, 'ylim', probYlims, 'ytick', {});
    title(['touch probability' labels{k}])
    xlabel('\itposition (m)')
    ylabel('\ittouch probability')


    % pimp out velocity figure
    subplot(2,3,2 + (k-1)*3)
    set(gca, 'xlim', [-obsPrePost(1) obsPrePost(2)], 'ylim', ylims)
    
    title(['velocity' labels{k}]);
    xlabel('\itposition (m)')
    ylabel('\itvelocity (m/s)')
    y1 = frameEdges(1)-obsPos;
    y2 = frameEdges(2)-obsPos;
    line([y1 y1], ylims, 'color', cmap{k}(1,:), 'linewidth', 2)
    line([y2 y2], ylims, 'color', cmap{k}(1,:), 'linewidth', 2)


    % pimp out bar graph
    subplot(2,3,3 + (k-1)*3)
    set(gca, 'ylim', [0 1], 'xlim', [.25 length(sessions)+.75])
    title(['success rate' labels{k}])
    xlabel('\itsession #')
    ylabel('\itfraction avoided')
end


% save fig
savefig(['obsAvoidance\figs\' mouse 'lights.fig'])
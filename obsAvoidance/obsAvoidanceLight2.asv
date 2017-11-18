function obsAvoidanceLight2(mouse, expName)

% compare obstacle avoidance with and without the obstacle light on
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
% trialRange = [0 .8]; % only include trials #s between these limits, so performance metrices are not dragged down when the animal is warming up or sated near the end
obsPos = .382; % m, position at which obstacle is in the middle of the frame // use getFrameTimes function to determine this value
frameEdges = [.336 .444]; % (m)
sig = .0025; % sigma for gaussian kernal
probYlims = [0 .004];
minTouchTime = 0; % only touches count that are >= minTouchTime
conditionLabels = {'light', 'no light'};


% initializations
sessionInfo = readtable([dataDir 'sessionInfo.xlsx']);

sessionInds = strcmp(sessionInfo.mouse, mouse) &...
              cellfun(@(x) any(strcmp(x, expName)), sessionInfo.experiment) &...
              sessionInfo.include;
sessions = sessionInfo.session(sessionInds);

posInterp = -obsPrePost(1) : posRes : obsPrePost(2);            % velocities will be interpolated across this grid of positional values
touchPosInterp = frameEdges(1): touchPosRes : frameEdges(2);    % positions for touch probability figs with be interpolated across this grid

gausKernel = arrayfun(@(x) (1/(sig*sqrt(2*pi))) * exp(-.5*(x/sig)^2), -sig*5:1/(1/touchPosRes):sig*5); % kernal for touch probability figure
gausKernel = gausKernel/sum(gausKernel);

data = struct(); % stores trial data for all sessions
dataInd = 0;

cmap = winter(length(sessions)*2);
cmaps{1} = cmap(end-length(sessions)+1:end,:);
cmaps{2} = cmap(1:length(sessions),:);




% COMPILE DATA

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
  
    
    % compute velocity
    vel = getVelocity(wheelPositions, .5, targetFs);
    
    
    % iterate over all trials
    for j = 1:length(obsOnTimes)
        
        dataInd = dataInd + 1;
        data(dataInd).session = sessions{1};
        data(dataInd).sessionNum = i;
        data(dataInd).name = mouse;
        
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
        data(dataInd).velocity = trialVelInterp;
        
        % find whether and where obstacle was toucheed
        trialTouchInds = touchOnTimes>obsOnTimes(j) & touchOnTimes<obsOffTimes(j);
        data(dataInd).avoided = ~any(trialTouchInds);
        trialTouchPositions = touchPositions(trialTouchInds);
        
        % compute touch probability as function of position
        touchCounts = histcounts(trialTouchPositions, touchPosInterp);
        touchCounts = touchCounts / length(trialTouchPositions);
        trialTouchProb = conv(touchCounts, gausKernel, 'same');
        trialTouchProb(isnan(trialTouchProb)) = 0;
        data(dataInd).touchProbability = trialTouchProb;
        
        % find whether light was on
        data(dataInd).obsLightOn = min(abs(obsOnTimes(j) - obsLightOnTimes)) < .5;
        
        % record position at which obstacle turned on
        data(dataInd).obsOnPositions = obsPos - obsOnPos;
    end
end



% PLOT EVERYTHING

% prepare figure
figure('name', mouse);


% plot touch probability
subplot(1,3,1)
allTouchProbs = reshape([data(:).touchProbability], length(data(1).touchProbability), length(data))';
touchPosCenters = touchPosInterp(1:end-1) - .5*touchPosRes - obsPos;
plot(touchPosCenters, mean(allTouchProbs([data.obsLightOn],:), 1), 'color', mean(cmaps{1},1), 'linewidth', 3); hold on;
plot(touchPosCenters, mean(allTouchProbs(~[data.obsLightOn],:), 1), 'color', mean(cmaps{2},1), 'linewidth', 3);

title('touch probability')
xlabel('\itposition (m)')
ylabel('\ittouch probability')

set(gca, 'xdir', 'reverse', 'xlim', [touchPosCenters(1) touchPosCenters(end)], 'ylim', probYlims, 'ytick', {});

legend(conditionLabels); legend('boxoff');



% plot velocity

subplot(1,3,2)
allVelocities = reshape([data(:).velocity], length(data(1).velocity), length(data))';

for i = 1:length(sessions)
    sessionInds = [data.sessionNum]==i;
    
    plot(posInterp, nanmean(allVelocities([data.obsLightOn] & sessionInds,:), 1), 'color', cmaps{1}(i,:), 'linewidth', 1); hold on
    plot(posInterp, nanmean(allVelocities(~[data.obsLightOn] & sessionInds,:), 1), 'color', cmaps{2}(i,:), 'linewidth', 1);
end

plot(posInterp, nanmean(allVelocities([data.obsLightOn],:), 1), 'color', mean(cmaps{1},1), 'linewidth', 3); hold on
plot(posInterp, nanmean(allVelocities(~[data.obsLightOn],:), 1), 'color', mean(cmaps{2},1), 'linewidth', 3);

title('velocity');
xlabel('\itposition (m)')
ylabel('\itvelocity (m/s)')

set(gca, 'xlim', [-obsPrePost(1) obsPrePost(2)], 'ylim', ylims)
x1 = frameEdges(1)-obsPos;
x2 = frameEdges(2)-obsPos;
line([x1 x1], ylims, 'color', [0 0 0], 'linewidth', 1)
line([x2 x2], ylims, 'color', [0 0 0], 'linewidth', 1)
obsOnPos = -mean([data.obsOnPositions]);
line([obsOnPos obsOnPos], ylims, 'color', [0 0 0], 'linewidth', 1)


% plot obstacle avoidance
subplot(1,3,3)

% compute avoidance per session for light on and off conditions
lightOnAvoidance  = nan(1,length(sessions));
lightOffAvoidance = nan(1,length(sessions));

for i = 1:length(sessions)
        
    trialOn = [data.sessionNum]==i & [data.obsLightOn];
    trialOff = [data.sessionNum]==i & ~[data.obsLightOn];
    
    lightOnAvoidance(i)  = sum([data(trialOn).avoided]) / sum(trialOn);
    lightOffAvoidance(i) = sum([data(trialOff).avoided]) / sum(trialOff);
end

scatter(ones(1,length(sessions))*1.5, lightOnAvoidance, 75, cmaps{1}, 'filled', 'jitter', 'on'); hold on
scatter(ones(1,length(sessions))*3.5, lightOffAvoidance, 75, cmaps{2}, 'filled', 'jitter', 'on');
line([1 2], repmat(mean(lightOnAvoidance),1,2), 'color', 'black', 'linewidth', 2);
line([3 4], repmat(mean(lightOffAvoidance),1,2), 'color', 'black', 'linewidth', 2);
set(gca, 'ylim', [0 1], 'xlim', [.5 4.5], 'xtick', [1.5 3.5], 'xticklabel', conditionLabels);

title('success rate')
xlabel('\itcondition')
ylabel('\itfraction avoided')

pimpFig;
set(gcf, 'menubar', 'none',...
         'units', 'inches',...
         'position', [4 4 14.5 4.5]);



% save fig
savefig(['obsAvoidance\figs\' mouse 'lights.fig'])

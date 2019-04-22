function plotSingleSessionVel(session)

% this figure is intended to show the trial strucutre by plotting vel as
% fcn of position relative to reward // position of 3 osbtacle will be
% overlaid on plot // single trials with be plotted over session mean

% TO DO: ***change line to wisk contact position, not nose position

% settings
% session = '180714_004';
posLims = [-5.4 0]; % (m) how much to plot before and after
posTics = 200; % how many elements in x axis position vector
speedTime = .01; % compute velocity over this interval
trialsToShow = 15; % how many individual trials to plot behind the mean
errorFcn = @(x) nanstd(x);
obsOnAlpha = .1;
trialAlpha = .4;
smoothing = .5; % (m)
topArea = .2; % add area to plot above to put lines indicating water reward and obstacle periods, with text // expressed as fraction of y range

obsOnColor = [0 0 0];
meanColor = [48 135 227]*.75 / 255;
% trialColors = lines(trialsToShow);
trialColors = repmat([1 1 1]*.4, trialsToShow,1);



% initializations
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'wheelPositions', 'wheelTimes', 'obsOnTimes', 'obsOffTimes', 'obsPositionsFixed', 'obsTimes', 'rewardTimes');
xPosits = linspace(posLims(1), posLims(2), posTics);
wheelVel = getVelocity(wheelPositions, speedTime, 1000);
rewardTimes = rewardTimes(2:end); % remove first reward, because can't trace backwards from it
smoothSmps = round(smoothing / mean(diff(xPosits)));

vels = nan(length(rewardTimes), length(xPosits));
obsOnPositions = interp1(wheelTimes, wheelPositions, obsOnTimes);
obsOffPositions = interp1(wheelTimes, wheelPositions, obsOffTimes);


trialObsOnPosits = nan(length(rewardTimes), 3); % 3 obstacles per reward
trialObsOffPosits = nan(length(rewardTimes), 3);
atNosePosits = nan(length(rewardTimes),3);

% collect vel vs position for all trials
for i = 1:length(rewardTimes)
    
    % get trial positions and velocity
    rewardPos = wheelPositions(find(wheelTimes>=rewardTimes(i),1,'first'));
    
    % get obstacle on, off, and at nose positions
    obsOnInds = find(obsOnPositions<rewardPos & obsOnPositions>rewardPos+posLims,3,'last');
    obsOffInds = find(obsOnPositions<rewardPos & obsOnPositions>rewardPos+posLims,3,'last');
    
    if length(obsOnInds==3) && isequal(obsOnInds, obsOffInds)
        
        trialObsOnPosits(i,:) = obsOnPositions(obsOnInds) - rewardPos;
        trialObsOffPosits(i,:) = obsOffPositions(obsOffInds)-rewardPos;
        
        trialPos = wheelPositions - rewardPos; % shift positions s.t. 0 is reward position
        bins = trialPos>posLims(1) & trialPos<posLims(2);
        trialPos = trialPos(bins);
        trialVel = wheelVel(bins);

        % remove duplicate positional values (would be better to average all values in a particular bin)
        [trialPos, uniqueInds] = unique(trialPos, 'stable');
        trialVel = trialVel(uniqueInds);
        
        % interpolate velocities across positional grid and store results
        vels(i,:) = smooth(interp1(trialPos, trialVel, xPosits, 'linear'), smoothSmps);
        
        % find positions at which obstacle reaches nose
        for j = 1:3
            atNoseTime = obsTimes(find(obsPositionsFixed>=0 & obsTimes>obsOnTimes(obsOnInds(j)),1,'first'));
            atNosePosits(i,j) = interp1(wheelTimes, wheelPositions, atNoseTime)-rewardPos;
        end
    end
end


% make plot
close all
figure('menubar', 'none', 'color', 'white', 'position', [2000 100 1200 300]); 

% plot trials
trials = randperm(length(rewardTimes), trialsToShow);
for i = 1:length(trials)
    plot(xPosits, vels(trials(i),:), 'color', [trialColors(i,:) trialAlpha]); hold on
end
yLims = max(0, get(gca, 'YLim'));
vertOffset = +.2*topArea*range(yLims); % offset of schematic lines above figure // change first term in expression to make it higher or lower

% add shaded areas where obstacle turns on
for i = 1:3
    
    % shaded area for obs on
    x = [nanmean(trialObsOnPosits(:,i),1) nanmean(trialObsOffPosits(:,i),1)];
    rectangle('Position', [x(1) yLims(1) diff(x) diff(yLims)], ...
        'FaceColor', [obsOnColor obsOnAlpha], 'EdgeColor', 'none');
    
    % line above obs on
    line([x(1) x(2)], [yLims(2) yLims(2)]+vertOffset, 'color', obsOnColor, 'linewidth', 3);
    if i==1; text(x(1), yLims(2)+2*vertOffset, 'obstacle enngaged', 'VerticalAlignment', 'bottom'); end
%     if i==1; text(x(2), yLims(2)+.5*topArea*range(yLims), 'obstacle enngaged'); end
    
    % line for obs at nose
    x = nanmean(atNosePosits(:,i),1);
    line([x x], yLims, 'color', obsOnColor, 'linewidth', 2);
end

% add water reward
line([xPosits(end) xPosits(end)], [yLims(2)+vertOffset yLims(2)+3*vertOffset], 'color', [0 0 .8], 'linewidth', 3)
text(xPosits(end), yLims(2)+4*vertOffset, 'water', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
% keyboard

% plot mean
shadedErrorBar(xPosits, vels, {@nanmean, errorFcn}, ...
        'lineprops', {'linewidth', 3, 'color', meanColor}); hold on;
% plot(xPosits, nanmean(vels,1), 'linewidth', 3, 'color', meanColor)




    
% pimp fig
set(gca, 'box', 'off', 'XLim', [posLims(1) posLims(2)], 'YLim', [yLims(1) yLims(2)+range(yLims)*topArea], 'TickDir', 'out', ...
    'XTick', -fliplr(0:1:-posLims(1)), 'YTick', 0:.2:yLims(2))
xlabel('distance to water reward (m)')
ylabel('velocity (m/s)')
line([posLims(1) posLims(1)], [yLims(2) yLims(2)+range(yLims)*topArea], 'color', 'white', 'linewidth', 3) % cover the extra vertical portion of the y axis with a white line - yes, this is a hack



% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'singleSessionVel');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');

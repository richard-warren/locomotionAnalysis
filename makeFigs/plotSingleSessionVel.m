function plotSingleSessionVel(session, varargin)

% this figure is intended to show the trial strucutre by plotting vel as
% fcn of position relative to reward // position of 3 osbtacle will be
% overlaid on plot // single trials with be plotted over session mean

% TO DO: ***change line to wisk contact position, not nose position

% settings
s.posLims = [-5.4 0]; % (m) how much to plot before and after
s.posTics = 200; % how many elements in x axis position vector
s.speedTime = .01; % compute velocity over this interval
s.trialsToShow = 15; % how many individual trials to plot behind the mean
s.errorFcn = @(x) nanstd(x);
s.obsOnAlpha = .1;
s.trialAlpha = .4;
s.smoothing = .5; % (m)
s.topArea = .2; % add area to plot above to put lines indicating water reward and obstacle periods, with text // expressed as fraction of y range
s.waterColor = 'blue';  % color of tick mark showing time of water delivery

s.obsOnColor = [0 0 0];
s.meanColor = [48 135 227]*.75 / 255;
s.trialColors = 'jet';



% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'wheelPositions', 'wheelTimes', 'obsOnTimes', 'obsOffTimes', 'obsPositionsFixed', 'obsTimes', 'rewardTimes');
xPosits = linspace(s.posLims(1), s.posLims(2), s.posTics);
wheelVel = getVelocity(wheelPositions, s.speedTime, 1000);
rewardTimes = rewardTimes(2:end); % remove first reward, because can't trace backwards from it
smoothSmps = round(s.smoothing / mean(diff(xPosits)));

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
    obsOnInds = find(obsOnPositions<rewardPos & obsOnPositions>rewardPos+s.posLims,3,'last');
    obsOffInds = find(obsOnPositions<rewardPos & obsOnPositions>rewardPos+s.posLims,3,'last');
    
    if length(obsOnInds)==3 && isequal(obsOnInds, obsOffInds)
        trialObsOnPosits(i,:) = obsOnPositions(obsOnInds) - rewardPos;
        trialObsOffPosits(i,:) = obsOffPositions(obsOffInds)-rewardPos;
        
        trialPos = wheelPositions - rewardPos; % shift positions s.t. 0 is reward position
        bins = trialPos>s.posLims(1) & trialPos<s.posLims(2);
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
figure('menubar', 'none', 'color', 'white', 'position', [2000 50 600 300]); 

% plot trials
trials = randperm(length(rewardTimes), s.trialsToShow);
for i = 1:length(trials)
    plot(xPosits, vels(trials(i),:), 'color', [s.trialColors(i,:) s.trialAlpha]); hold on
end
yLims = max(0, get(gca, 'YLim'));
vertOffset = +.2*s.topArea*range(yLims); % offset of schematic lines above figure // change first term in expression to make it higher or lower

% add shaded areas where obstacle turns on
for i = 1:3
    
    % shaded area for obs on
    x = [nanmean(trialObsOnPosits(:,i),1) nanmean(trialObsOffPosits(:,i),1)];
    rectangle('Position', [x(1) yLims(1) diff(x) diff(yLims)], ...
        'FaceColor', [s.obsOnColor s.obsOnAlpha], 'EdgeColor', 'none');
    
    % line above obs on
    line([x(1) x(2)], [yLims(2) yLims(2)]+vertOffset, 'color', s.obsOnColor, 'linewidth', 3);
    if i==1; text(x(1), yLims(2)+2*vertOffset, 'obstacle engaged', 'VerticalAlignment', 'bottom'); end
    
    % line for obs at nose
    x = nanmean(atNosePosits(:,i),1);
    line([x x], yLims, 'color', get(gca, 'xcolor'), 'linewidth', 2);
end

% add water reward
line([xPosits(end) xPosits(end)], [yLims(2)+vertOffset yLims(2)+3*vertOffset], 'color', s.waterColor, 'linewidth', 3)
text(xPosits(end), yLims(2)+4*vertOffset, 'water', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')

% plot mean
shadedErrorBar(xPosits, vels, {@nanmean, s.errorFcn}, ...
        'lineprops', {'linewidth', 3, 'color', s.meanColor}); hold on;




    
% pimp fig
set(gca, 'box', 'off', 'XLim', [s.posLims(1) s.posLims(2)], 'YLim', [yLims(1) yLims(2)+range(yLims)*s.topArea], 'TickDir', 'out', ...
    'XTick', -fliplr(0:1:-s.posLims(1)), 'YTick', 0:.2:yLims(2))
xlabel('distance to water reward (m)')
ylabel('velocity (m/s)')
line([s.posLims(1) s.posLims(1)], [yLims(2) yLims(2)+range(yLims)*s.topArea], 'color', 'white', 'linewidth', 3) % cover the extra vertical portion of the y axis with a white line - yes, this is a hack

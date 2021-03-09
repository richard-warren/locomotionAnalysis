% heatmaps for all units showing responses to entire trial structure

unitInfo = getUnitInfo('nucleiOnly', true);
sessions = unique(unitInfo.session);
paper2_config;

%% compute median times of events for each session (all relative to preceding reward)
% events are: reward_pevious, (obsOn, wisk, obsOff)x3, reward

% settings
preObsOn = 2;    % (s) create initial event preObsOn s before first obs on
postReward = 5;  % (s) create final event postReward s after reward
tres = .01;      % (s) temporal resolution


inits = repelem({nan(length(sessions),1)}, 12);
events = table(inits{:}, 'RowNames', sessions', ...
    'VariableNames', {'start', 'on1', 'wisk1', 'off1', 'on2', 'wisk2', 'off2', 'on3', 'wisk3', 'off3', 'reward', 'end'});

for i = 1:length(sessions)
    
    % load session data
    d = load(['E:\lab_files\paper2\modelling\runAnalyzed\' sessions{i} '_runAnalyzed.mat'], ...
        'obsOnTimes', 'obsOffTimes', 'wiskContactTimes', 'rewardTimes');
    firstReward = max(2, find(d.rewardTimes>d.obsOnTimes(1), 1, 'first'));  % first reward with preceding obstacles (or second reward if first reward has preceding obstacles)
    times = nan(length(d.rewardTimes), 12);  % matrix of times of all events
    
%     if strcmp(sessions{i}, '201216_000'); keyboard; end
    
    for j = firstReward:length(d.rewardTimes)
        
        % get times of 3 obs on, 3 obs off, and 3 wisk contact times for trial
        obson = d.obsOnTimes(d.obsOnTimes>d.rewardTimes(j-1) & d.obsOnTimes<d.rewardTimes(j));
        obsoff = d.obsOffTimes(d.obsOffTimes>d.rewardTimes(j-1) & d.obsOffTimes<d.rewardTimes(j));
        wisk = d.wiskContactTimes(d.wiskContactTimes>d.rewardTimes(j-1) & d.wiskContactTimes<d.rewardTimes(j));
        
        if length(obson)==3 && length(obsoff)==3 && length(wisk)==3
            times(j, 1)  = obson(1) - preObsOn;
            times(j, 11) = d.rewardTimes(j);
            times(j, 12) = d.rewardTimes(j) + postReward;
            times(j, [2 5 8]) = obson;
            times(j, [3 6 9]) = wisk;
            times(j, [4 7 10]) = obsoff;
        end
    end
    times = times - times(:,1);
    events{sessions{i}, :} = nanmedian(times, 1);
end

meanEventTimes = mean(events{:,:}, 1);  % average of median times of 10 events across mice

%% for each session, compute warped time axis

% inter-event duractions are stretched to matched averages cmputed above,
% and each reward interval starts at 0, such that average firing rates can
% be computed easily wrt this time axis


inits = repelem({cell(length(sessions),1)}, 6);
data = cat(2, events, table(inits{:}, ...
    'VariableNames', {'vel', 'lickRate', 'spkRates', 'unit_ids', 't', 't_stretched'}));


tic
for i = 1:length(sessions)
    disp(i/length(sessions))
    
    % load session data
    d = load(['E:\lab_files\paper2\modelling\runAnalyzed\' sessions{i} '_runAnalyzed.mat'], ...
        'obsOnTimes', 'obsOffTimes', 'wiskContactTimes', 'rewardTimes', ...
        'wheelPositions', 'wheelTimes', 'lickTimes', 'isLightOn');
    e = load(['E:\lab_files\paper2\modelling\neuralData\' sessions{i} '_neuralData.mat']);
    t = d.wheelTimes(1) : tres : d.wheelTimes(end);
    
    vel = getVelocity(d.wheelPositions, .05, 1/nanmedian(diff(d.wheelTimes)));
    vel = interp1(d.wheelTimes, vel, t);
    lickRate = getFiringRate(d.lickTimes, 'kernel', 'gauss', 'kernelSig', .05, 'times', t);
    spkRates = interp1(e.timeStamps, e.spkRates', t)'; % put fr on same time axis as lick rate and vel
    
%     % temp (mask light on trials)
%     epochs = d.obsOnTimes(d.isLightOn) + [0 2];
%     lightOnBins = any(t>epochs(:,1) & t<epochs(:,2), 1);
%     vel(lightOnBins) = nan;
    
    t_stretched = nan(size(t));
    firstReward = max(2, find(d.rewardTimes>d.obsOnTimes(1), 1, 'first'));  % first reward with preceding obstacles (or second reward if first reward has preceding obstacles)
    for j = firstReward:length(d.rewardTimes)
        
        % get times of 3 obs on, 3 obs off, and 3 wisk contact times for trial
        obson = d.obsOnTimes(d.obsOnTimes>d.rewardTimes(j-1) & d.obsOnTimes<d.rewardTimes(j));
        obsoff = d.obsOffTimes(d.obsOffTimes>d.rewardTimes(j-1) & d.obsOffTimes<d.rewardTimes(j));
        wisk = d.wiskContactTimes(d.wiskContactTimes>d.rewardTimes(j-1) & d.wiskContactTimes<d.rewardTimes(j));
        
        times = nan(11,1);
        if length(obson)==3 && length(obsoff)==3 && length(wisk)==3
            times(1)  = obson(1) - preObsOn;
            times(11) = d.rewardTimes(j);
            times(12) = d.rewardTimes(j) + postReward;
            times([2 5 8]) = obson;
            times([3 6 9]) = wisk;
            times([4 7 10]) = obsoff;
            
            bins = t>=times(1) & t<=times(end);
            t_stretched(bins) = interp1(times, meanEventTimes, t(bins), 'linear');
        end
    end
    
    data.vel{i} = vel;
    data.lickRate{i} = lickRate;
    data.spkRates{i} = spkRates;
    data.unit_ids{i} = e.unit_ids;
    data.t{i} = t;
    data.t_stretched{i} = t_stretched;
end
toc

%% compute average vel, lick rate, firing rates per session

% settings
pts = 20;       % points along x axis
windowSz = .01;  % fraction of x axis


x = linspace(meanEventTimes(1), meanEventTimes(end), 1000);
data.vel_mean = nan(length(sessions), length(x));
data.lickrate_mean = nan(length(sessions), length(x));
unitInfo.frMean = nan(height(unitInfo), length(x));
win = range(x) * windowSz * [-.5 .5];  % sliding widndow range

for i = 1:length(sessions)
    disp(i/length(sessions))
    t = data{sessions{i}, 't_stretched'}{1};
    stacked = [data{sessions{i}, 'spkRates'}{1}; ...
               data{sessions{i}, 'vel'}{1}; ...
               data{sessions{i}, 'lickRate'}{1}];
    means = nan(size(stacked,1), length(x));
    
    for j = 1:length(x)
        bins = t>(x(j)+win(1)) & t<=(x(j)+win(2));
        if any(bins)
            means(:,j) = nanmean(stacked(:, bins), 2);
        end
    end
    
    data.vel_mean(i,:) = means(end-1,:);
    data.lickrate_mean(i,:) = means(end,:);
    
    unitInds = find(strcmp(unitInfo.session, sessions{i}) & ...
        ismember(unitInfo.unit, data{sessions{i}, 'unit_ids'}{1}));
    for j = 1:length(unitInds)
        sesBin = data{sessions{i}, 'unit_ids'}{1} == unitInfo.unit(unitInds(j));
        unitInfo.frMean(unitInds(j), :) = means(sesBin, :);
    end
end


%%

% settings
velLims = [0 .6];
lickLims = [0 10];
nclusters = 12;
pcs = 6;

close all
figure('position', [98.00 69.00 584.00 1216.00], 'color', 'white', 'menubar', 'none')
subplot(4,1,1); hold on


% events
wiskt = meanEventTimes(contains(events.Properties.VariableNames, 'wisk'));
plot(repelem(wiskt,2,1), velLims, '-', 'color', cfg.wiskColor, 'LineWidth', 2)  % wisk

rewardt = meanEventTimes(contains(events.Properties.VariableNames, 'reward'));
plot(repelem(rewardt,2,1), velLims, '-', 'color', cfg.lickColor, 'LineWidth', 2)  % wisk

ont  = meanEventTimes(contains(events.Properties.VariableNames, 'on'));
offt = meanEventTimes(contains(events.Properties.VariableNames, 'off'));
X = [ont; offt; offt; ont];
Y = repmat([velLims(2) velLims(2) velLims(1) velLims(1)]', 1, 3);
patch(X, Y, cfg.obsColor, 'EdgeColor', 'none', 'FaceAlpha', .1);

% vel
yyaxis left;
mn = mean(data.vel_mean,1);
st = std(data.vel_mean,[],1);
plot(x, mn, 'LineWidth', 2, 'Color', cfg.velColor);
patch([x fliplr(x)], [mn+st fliplr(mn-st)], cfg.velColor, 'EdgeColor', 'none', 'FaceAlpha', .1)
ylabel('velocity (m/s)')
set(gca, 'YColor', cfg.velColor, 'YLim', velLims, cfg.axArgs{:}); limitticks

% licks
yyaxis right;
mn = mean(data.lickrate_mean,1);
st = std(data.lickrate_mean,[],1);
plot(x, mn, 'LineWidth', 2, 'Color', cfg.lickColor);
patch([x fliplr(x)], [mn+st fliplr(mn-st)], cfg.lickColor, 'EdgeColor', 'none', 'FaceAlpha', .1)
ylabel('lick rate (licks/s)')
set(gca, 'XLim', [x(1) x(end)], 'YLim', lickLims, 'XColor', 'none', 'YColor', cfg.lickColor, cfg.axArgs{:}); limitticks

% heatmaps
lims = [-6 6];  % (z-score)

% zscore, remove outliers
resp = (unitInfo.frMean - nanmean(unitInfo.frMean,2)) ./ nanstd(unitInfo.frMean,[],2);
resp = resp(~any(isnan(resp),2),:);  % rows with no nans
resp(resp<lims(1)) = lims(1); resp(resp>lims(2)) = lims(2);
nrows = size(resp,1);

% cluster and sort by response type
groups = clusterResponses(resp, 'plot', false, 'nclusters', nclusters, 'pcs', pcs);
[~, maxind] = max(resp,[],2);

[~, sortInds] = sortrows([groups maxind], 'descend');

% plot
subplot(4,1,2:4); hold on
imagesc(x, 1:nrows, resp(sortInds,:), lims)
colormap(cfg.heatmapColors)

plot(repelem(wiskt,2,1), [.5 nrows+.5], '-', 'color', cfg.wiskColor, 'LineWidth', 2)  % wisk
plot(repelem(rewardt,2,1), [.5 nrows+.5], '-', 'color', cfg.lickColor, 'LineWidth', 2)  % wisk
Y = repmat([nrows+.5 nrows+.5 .5 .5]', 1, 3);
patch(X, Y, cfg.obsColor, 'EdgeColor', 'none', 'FaceAlpha', .1);

set(gca, 'XLim', [x(1) x(end)], 'YLim', [.5 nrows+.5], 'XColor', 'none', 'YColor', 'none', cfg.axArgs{:})

% add labels
text(wiskt(1), nrows, 'obstacle', 'Color', cfg.obsColor, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
text(rewardt(1), nrows, 'water reward', 'Color', cfg.lickColor, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
plot([0 2], [.5 .5], 'color', 'black', 'LineWidth', 2)
text(0, 0, '2 seconds', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')

saveas(gcf, 'E:\lab_files\paper2\paper_figures\matlab_exports\overallheatmaps', 'svg')


















%%
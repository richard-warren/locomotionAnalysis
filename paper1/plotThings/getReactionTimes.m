%% inits
fprintf('loading... '); load(fullfile(getenv('SSD'), 'paper1', 'baseline_data.mat'), 'data'); disp('baseline data loaded!')

% global settings
paper1_config;  % initialize global settings

flat = flattenData(data, ...
    {'mouse', 'session', 'trial', 'preModPawKin', 'modPawKin', 'isBigStep', 'isLightOn', 'velAtWiskContact', 'contactInd'});
mice = unique({flat.mouse});

% settings
neighborNum = 40; % how many similar trials to use to construct control distribution
dt = .004; % 1/fps of camera
xlims = [-10 80]; % (ms)
ylims = [-20 20];
minVel = .4;
modPawFractions = [.25 .75];  % only include trials where whisker contact occurs within these fraction limits of the total swing duration
colors = decisionColors;

%% initializations

% find where the swing ends in each step
lastStepInd = nan(1,length(flat));
for i = find(cellfun(@(x) ~all(isnan(x(:))), {flat.modPawKin})) % only search inds where modPawKin are computed
    lastStepInd(i) = find(any(diff(flat(i).modPawKin,[],2),1),1,'last');
end

% include trials in which whisker contact occurs mid-swing, light is off, and running is fast
midSwingBins = [flat.contactInd] > modPawFractions(1)*lastStepInd & [flat.contactInd] < modPawFractions(2)*lastStepInd;
flat = flat(midSwingBins & ...
            ~[flat.isLightOn] & ...
            [flat.velAtWiskContact]>minVel); % add conditionals here

swingMaxSmps = size(flat(1).modPawKin,2);  % inherits max samples for swing locations from getKinematicData
times = linspace(-swingMaxSmps*dt, swingMaxSmps*dt, swingMaxSmps*2+1);  % time relative to whisker contact // everything is centered around the moment of contact

% get kinematics
modLocationsRaw = permute(cat(3,flat.modPawKin),[3,1,2]);  % mod paw
controlLocationsRaw = cellfun(@(x) squeeze(x(end,:,:)), {flat.preModPawKin}, 'UniformOutput', false);  % step before mod paw
controlLocationsRaw = permute(cat(3,controlLocationsRaw{:}), [3,1,2]);
control2LocationsRaw = cellfun(@(x) squeeze(x(end-1,:,:)), {flat.preModPawKin}, 'UniformOutput', false);  % step before that
control2LocationsRaw = permute(cat(3,control2LocationsRaw{:}), [3,1,2]);

% subtract x, y, and z values s.t. they start at zero at the start of every trial
control2LocationsRaw = control2LocationsRaw - control2LocationsRaw(:,:,1);
controlLocationsRaw = controlLocationsRaw - controlLocationsRaw(:,:,1);
modLocationsRaw = modLocationsRaw - modLocationsRaw(:,:,1);



% inits
% these mats will contain kinematics centered at whisker contact
modLocations = nan(length(flat), 3, length(times));
controlLocations = nan(length(flat), 3, length(times), neighborNum);   % every mod step will have neighborNum matched control steps
control2Locations = nan(length(flat), 3, length(times), neighborNum);  % every best matched control step will ahve neighborNum matched control steps from the previous step
modDifs = nan(length(flat), 3, length(times));       % diffs between mod and control step
controlDifs = nan(length(flat), 3, length(times));   % diffs between control step and previous control steps



for i = 1:length(flat)
    disp(i/length(flat))
    
    % get temporal inds
    trialInds = (swingMaxSmps+1 : swingMaxSmps*2) - (flat(i).contactInd+1);  % everything is centered at the momemnt of contact, s.t. moment of contact occurs at time 0
    tbins = 1:flat(i).contactInd;  % inds of samples occuring before whisker contact
    
    % get nearest neighbor control trials
    % (match based on kinematics prior to whisker contact)
    
%     controlLocationsSub = squeeze(controlLocationsRaw(:,dimToMatch,1:flat(i).contactInd));
%     modLocationsSub = squeeze(modLocationsRaw(i,dimToMatch,1:flat(i).contactInd))';
    controlLocationsSub = [squeeze(controlLocationsRaw(:,1,tbins)) squeeze(controlLocationsRaw(:,3,tbins))];  % contact x and z
    modLocationsSub     = [squeeze(modLocationsRaw(i,1,tbins))'    squeeze(modLocationsRaw(i,3,tbins))'];  % contact x and z
    controlInds = knnsearch(controlLocationsSub, modLocationsSub, 'k', neighborNum); % inds of nearest trials in terms of velocity
    controlLocations(i,:,trialInds,:) = permute(controlLocationsRaw(controlInds,:,:), [2,3,1]);

    % get nearest neighbor control trials for best matched control trial
%     control2LocationsSub = squeeze(control2LocationsRaw(:,dimToMatch,1:flat(i).contactInd));
%     controlLocationsSub = squeeze(controlLocationsRaw(controlInds(1),dimToMatch,1:flat(i).contactInd))';
    control2LocationsSub = [squeeze(control2LocationsRaw(:,1,tbins))              squeeze(control2LocationsRaw(:,3,tbins))];  % contact x and z
    controlLocationsSub  = [squeeze(controlLocationsRaw(controlInds(1),1,tbins))' squeeze(controlLocationsRaw(controlInds(1),3,tbins))'];  % contact x and z
    control2Inds = knnsearch(control2LocationsSub, controlLocationsSub, 'k', neighborNum); % inds of nearest trials in terms of velocity
    control2Locations(i,:,trialInds,:) = permute(control2LocationsRaw(control2Inds,:,:), [2,3,1]);

    % get dif between each control trial and average of previous control trials
    control2Mean = nanmean(control2LocationsRaw(control2Inds,:,:), 1);
    controlDifs(i,:,trialInds) = controlLocations(i,:,trialInds) - control2Mean;

    % get mod locations
    modLocations(i,:,trialInds) = modLocationsRaw(i,:,:);

    % get dif between average of control locations and mod locations
    controlMean = nanmean(controlLocationsRaw(controlInds,:,:), 1);
    modDifs(i,:,trialInds) = modLocations(i,:,trialInds) - controlMean;
end



%% main plot

figure('color', 'white', 'menubar', 'none', 'position', [2.00 1097.00 1278.00 313.00]);

errorFcn = @(x) nanstd(x)/sqrt(size(x,1));
% errorFcn = @(x) nanstd(x);
dims = {'x', 'y', 'z'};

mouseDifs = nan(length(mice), 3, size(modDifs,3), 2); % [mouse X dim X kinDiffs X (big step, little step, control)]

for i = 1:length(mice)
    bins = strcmp({flat.mouse}, mice{i});
    mouseDifs(i,:,:,1) = nanmean(modDifs(bins & ~[flat.isBigStep],:,:),1);
    mouseDifs(i,:,:,2) = nanmean(modDifs(bins & [flat.isBigStep],:,:),1);
    mouseDifs(i,:,:,3) = nanmean(controlDifs(bins,:,:),1);
end

for i = 1:length(dims)
    
    % get difs for one dimension
    subplot(1,length(dims),i)    
    shadedErrorBar(times*1000, squeeze(mouseDifs(:,i,:,1)*1000), {@nanmean, errorFcn}, ...
        'lineprops', {'linewidth', 3, 'color', colors(1,:)}); hold on
    shadedErrorBar(times*1000, squeeze(mouseDifs(:,i,:,2)*1000), {@nanmean, errorFcn}, ...
        'lineprops', {'linewidth', 3, 'color', colors(2,:)});
    shadedErrorBar(times*1000, squeeze(mouseDifs(:,i,:,3)*1000), {@nanmean, errorFcn}, ...
        'lineprops', {'linewidth', 3, 'color', [.2 .2 .2]});

    % pimp fig
    set(gca, 'xlim', xlims, 'ylim', ylims)
    ylabel(['\Delta' dims{i} ' (mm)'])
    xlabel('time from whisker contact (ms)')
    line([0 0], get(gca,'ylim'), 'color', [0 0 0])
end

saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'reactionTimes.svg'), 'svg');




%% plot single trial dif from control population

close all; figure('color', 'white', 'menubar', 'none', 'position', [2.00 2.00 1278.00 1408.00]);

rows = 8;
cols = 4;
trialInds = randperm(length(flat), (rows*cols)/2);
xlimsTemp = [-.1 .2];
dim = 1;

errorFcn = @(x) 2*nanstd(x);

plotInd = 1;
for i = trialInds
    
    subaxis(rows, cols, plotInd, 'spacing', .02, 'margin', .08, 'padding', 0)
    if flat(i).isBigStep; colorInd=2; else; colorInd=1; end
    plot(times, squeeze(modLocations(i,dim,:)), 'color', colors(colorInd,:), 'linewidth', 3); hold on;
    shadedErrorBar(times, squeeze(controlLocations(i,dim,:,:))', {@nanmean, errorFcn}, ...
        'lineprops', {'linewidth', 3, 'color', [.5 .5 .5]});
    set(gca, 'box', 'off', 'xtick', times(1):.08:times(end), 'ytick', [], 'ycolor', 'white', 'xlim', xlimsTemp);
    line([0 0], get(gca,'ylim'), 'color', [0 0 0])
    plotInd = plotInd + 1;
    
    subaxis(rows, cols, plotInd, 'spacing', .02, 'margin', .08, 'padding', 0)
    plot(times, squeeze(controlLocations(i,dim,:,1)), 'color', [0 0 0], 'linewidth', 3); hold on;
    shadedErrorBar(times, squeeze(control2Locations(i,dim,:,:))', {@nanmean, errorFcn}, ...
        'lineprops', {'linewidth', 3, 'color', [.5 .5 .5]});
    set(gca, 'box', 'off', 'xtick', times(1):.08:times(end), 'ytick', [], 'ycolor', 'white', 'xlim', xlimsTemp);
    plotInd = plotInd + 1;
    
    
    line([0 0], get(gca,'ylim'), 'color', [0 0 0])
    
end


% saveas(gcf, [getenv('OBSDATADIR') 'figures\reactionTimes\reactionTimesTrials.png']);



%% plot for each mouse separately


close all; figure('color', 'white', 'menubar', 'none', 'position', [2.00 2.00 1278.00 1408.00]);
errorFcn = @(x) nanstd(x)/sqrt(size(x,1));
% errorFcn = @(x) nanstd(x);
dims = {'x', 'y', 'z'};

for mouse = 1:length(mice)
    for i = 1:length(dims)
        
        subplot(ceil(length(mice)/2), 6, (mouse-1)*3 + i)
    
        % get diffs for one dimension
        bins = strcmp({flat.mouse}, mice{mouse});
        controlDifsSub = squeeze(controlDifs(bins,i,:));
        modDifsSub = squeeze(modDifs(bins,i,:));

        try
            shadedErrorBar(times, -controlDifsSub*1000, {@nanmean, errorFcn}, ...
                'lineprops', {'linewidth', 3, 'color', [0 0 0]}); hold on;
            shadedErrorBar(times, modDifsSub(~[flat(bins).isBigStep],:)*1000, {@nanmean, errorFcn}, ...
                'lineprops', {'linewidth', 3, 'color', colors(1,:)});
            shadedErrorBar(times, modDifsSub([flat(bins).isBigStep],:)*1000, {@nanmean, errorFcn}, ...
                'lineprops', {'linewidth', 3, 'color', colors(2,:)});
        end

        % pimp fig
        set(gca, 'xlim', xlimsTemp, 'ylim', ylims*2)
        ylabel(['\Delta' dims{i} ' (mm)'])
        xlabel('time (ms)')
        line([0 0], get(gca,'ylim'), 'color', [0 0 0])
    end
end




%% figure out the god damn latencies!
% (1) see when kinematics exceed percentile thresholds of control distribution

percentiles = [5 95];
stdThresh = 2.5;
dim = 3;
usePercentile = false;

zeroInd = find(times==0);
latencies = nan(1,length(flat));
latenciesControl = nan(1,length(flat));

for i = 1:length(flat)
    
    % mod
    modTrial = squeeze(modLocations(i,dim,:));
    
    % percentile lims
    if usePercentile
        lims = prctile(squeeze(controlLocations(i,dim,:,:))', percentiles, 1)';
    
    % std lims
    else
        ctlMean = nanmean(squeeze(controlLocations(i,dim,:,:))');
        ctlStd = std(squeeze(controlLocations(i,dim,:,:))');
        lims = ctlMean' + stdThresh*[-ctlStd; ctlStd]';
    end
    
    firstDeviatedInd = find((modTrial<lims(:,1) | modTrial>lims(:,2)) & times'>0, 1, 'first');
%     firstDeviatedInd = find(modTrial<lims(:,1) | modTrial>lims(:,2), 1, 'first');
    if ~isempty(firstDeviatedInd)
        latencies(i) = times(firstDeviatedInd);
    end
    
    
    % control
    controlTrial = squeeze(controlLocations(i,dim,:,1));  % control trial best matching modified trial (the first index in the final dimension)
    
    % percentile lims
    if usePercentile
        lims = prctile(squeeze(control2Locations(i,dim,:,:))', percentiles, 1)';
    
    % std lims
    else
        ctlMean = nanmean(squeeze(control2Locations(i,dim,:,:))');
        ctlStd = std(squeeze(control2Locations(i,dim,:,:))');
        lims = ctlMean' + stdThresh*[-ctlStd; ctlStd]';
    end
    
    firstDeviatedInd = find((controlTrial<lims(:,1) | controlTrial>lims(:,2)) & times'>0, 1, 'first');
%     firstDeviatedInd = find(controlTrial<lims(:,1) | controlTrial>lims(:,2), 1, 'first');
    if ~isempty(firstDeviatedInd)
        latenciesControl(i) = times(firstDeviatedInd);
    end
end

close all; figure('color', 'white', 'position', [550.00 715.00 560.00 420.00]);
subplot(2,1,1); histogram(latencies, 20); set(gca, 'xlim', [-.1 .1])
subplot(2,1,2); histogram(latenciesControl, 20); set(gca, 'xlim', [-.1 .1])

fprintf('latency estimate: %.2f ms\n', nanmedian(latencies)*1000)

% saveas(gcf, [getenv('OBSDATADIR') 'figures\reactionTimes\reactionTimeDistributions.png']);

%% (2) sliding window statistical tests

binEdges = -.1:.005:.1;  % s
dim = 1;
pvals = nan(1, length(binEdges)-1);

for i = 1:length(binEdges)-1
    tbins = times>=binEdges(i) & times<binEdges(i+1);
    ctl = nanmean(squeeze(controlDifs(:,dim,tbins)),2);  % avg control dif for each trial in this time bin
    mod = nanmean(squeeze(modDifs(:,dim,tbins)),2);
    [~, p] = ttest(ctl, mod);
    pvals(i) = p;
end

binCenters = binEdges(1:end-1) + diff(binEdges(1:2));
close all; figure('color', 'white', 'position', [550.00 715.00 560.00 420.00]);
scatter(binCenters, pvals, [], 'black', 'filled');
set(gca, 'ylim', [0 1], 'xlim', binEdges([1 end]))










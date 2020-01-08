% function getReactionTimes(data)

flat = getNestedStructFields(data, ...
    {'mouse', 'session', 'trial', 'preModPawKin', 'modPawKin', 'isBigStep', 'isLightOn', 'modPawContactInd', 'velAtWiskContact'});


% settings
neighborNum = 20; % how many similar trials to use to construct control distribution
dt = .004; % 1/fps of camera
xlims = [-10 80]; % (ms)
ylims = [-20 20];
colors = [.25 1 1; .25 1 .25];
dimToMatch = 3;
minVel = .4;
modPawLims = [5, 20]; % only include trials where contact occurs within modPawLims samples of swing start


% initializations

% find where the swing ends in each step
% lastStepInd = nan(1,length(flat));
% for i = find(cellfun(@(x) ~all(isnan(x(:))), {flat.modPawKin})) % only search inds where modPawKin are computed
%     lastStepInd(i) = find(any(diff(flat(i).modPawKin,[],2),1),1,'last');
% end

flat = flat([flat.modPawContactInd]>=modPawLims(1) & [flat.modPawContactInd]<=modPawLims(2) & ...
            ~[flat.isLightOn] & ...
            [flat.velAtWiskContact]>minVel); % add conditionals here

swingMaxSmps = size(flat(1).modPawKin,2); % inherits max samples for swing locations from getKinematicData
times = linspace(-swingMaxSmps*dt, swingMaxSmps*dt, swingMaxSmps*2+1); % everything is centered around the moment of contact

% get kinematics
modLocationsRaw = permute(cat(3,flat.modPawKin),[3,1,2]);
controlLocationsRaw = cellfun(@(x) squeeze(x(end,:,:)), {flat.preModPawKin}, 'UniformOutput', false);
controlLocationsRaw = permute(cat(3,controlLocationsRaw{:}), [3,1,2]);
control2LocationsRaw = cellfun(@(x) squeeze(x(end-1,:,:)), {flat.preModPawKin}, 'UniformOutput', false);
control2LocationsRaw = permute(cat(3,control2LocationsRaw{:}), [3,1,2]);


% subtract x and y values s.t. they start at zero at the start of every trial
control2LocationsRaw = control2LocationsRaw - repmat(control2LocationsRaw(:,:,1), 1, 1, swingMaxSmps);
controlLocationsRaw = controlLocationsRaw - repmat(controlLocationsRaw(:,:,1), 1, 1, swingMaxSmps);
modLocationsRaw = modLocationsRaw - repmat(modLocationsRaw(:,:,1), 1, 1, swingMaxSmps);



% initialization data containers
modLocations = nan(length(flat), 3, length(times));
controlLocations = nan(length(flat), 3, length(times), neighborNum);
control2Locations = nan(length(flat), 3, length(times), neighborNum);
modDifs = nan(length(flat), 3, length(times));
controlDifs = nan(length(flat), 3, length(times));


for i = 1:length(flat)
%     fprintf('\b\b\b\b\b')
%     fprintf('%.3f', i/length(flat))
    disp(i)
    
    % get temporal inds
    trialInds = (swingMaxSmps+1 : swingMaxSmps*2) - (flat(i).modPawContactInd+1); % everything is centered at the momemnt of contact, s.t. moment of contact occurs at time 0

    % get nearest neighbor control trials
    controlLocationsSub = squeeze(controlLocationsRaw(:,dimToMatch,1:flat(i).modPawContactInd));
    modLocationsSub = squeeze(modLocationsRaw(i,dimToMatch,1:flat(i).modPawContactInd))';
    controlInds = knnsearch(controlLocationsSub, modLocationsSub, 'k', neighborNum); % inds of nearest trials in terms of velocity
    controlLocations(i,:,trialInds,:) = permute(controlLocationsRaw(controlInds,:,:), [2,3,1]);

    % get nearest neighbor control trials for best control trial
    control2LocationsSub = squeeze(control2LocationsRaw(:,dimToMatch,1:flat(i).modPawContactInd));
    controlLocationsSub = squeeze(controlLocationsRaw(controlInds(1),dimToMatch,1:flat(i).modPawContactInd))';
    control2Inds = knnsearch(control2LocationsSub, controlLocationsSub, 'k', neighborNum); % inds of nearest trials in terms of velocity
    control2Locations(i,:,trialInds,:) = permute(control2LocationsRaw(control2Inds,:,:), [2,3,1]);

    % get dif between average of remaining control locations and control trial
    controlMean = nanmean(control2LocationsRaw(control2Inds,:,:), 1);
%     controlDifs(i,:,trialInds) = abs(controlLocations(i,:,trialInds) - controlMean);
    controlDifs(i,:,trialInds) = controlLocations(i,:,trialInds) - controlMean;

    % get mod locations and dif from control
    modLocations(i,:,trialInds) = modLocationsRaw(i,:,:);

    % get dif between average of control locations and mod locations
    controlMean = nanmean(controlLocationsRaw(controlInds,:,:), 1);
%     modDifs(i,:,trialInds) = abs(modLocations(i,:,trialInds) - controlMean);
    modDifs(i,:,trialInds) = modLocations(i,:,trialInds) - controlMean;
end



%% old version of the plot

figure('color', 'white', 'menubar', 'none', 'position', [2092 554 1631 400], 'InvertHardcopy', 'off');

errorFcn = @(x) nanstd(x)/sqrt(size(x,1));
% errorFcn = @(x) nanstd(x);
dims = {'x', 'y', 'z'};

for i = 1:length(dims)
    
    % get diifs for one dimension
    controlDifsSub = squeeze(controlDifs(:,i,:));
    modDifsSub = squeeze(modDifs(:,i,:));
    
    
    subplot(1,length(dims),i)
    
    shadedErrorBar(times*1000, controlDifsSub*1000, {@nanmean, errorFcn}, ...
        'lineprops', {'linewidth', 3, 'color', [0 0 0]}); hold on;
%     shadedErrorBar(times*1000, -controlDifsSub*1000, {@nanmean, errorFcn}, ...
%         'lineprops', {'linewidth', 3, 'color', [0 0 0]}); hold on;
    
    shadedErrorBar(times*1000, modDifsSub(~[flat.isBigStep],:)*1000, {@nanmean, errorFcn}, ...
        'lineprops', {'linewidth', 3, 'color', colors(1,:)});
    shadedErrorBar(times*1000, modDifsSub([flat.isBigStep],:)*1000, {@nanmean, errorFcn}, ...
        'lineprops', {'linewidth', 3, 'color', colors(2,:)});

    % pimp fig
    set(gca, 'xlim', xlims, 'ylim', ylims)
    ylabel(['\Delta' dims{i} ' (mm)'])
    xlabel('time (ms)')
    line([0 0], get(gca,'ylim'), 'color', [0 0 0])
end

blackenFig
% saveas(gcf, [getenv('OBSDATADIR') 'figures\reactionTimes\reactionTimes.png']);
% savefig(fullfile(getenv('OBSDATADIR'), 'figures\reactionTimes\reactionTimes.fig'))

%% subtracting control diffs, then measuring variability across mice!

figure('color', 'white', 'menubar', 'none', 'position', [2092 554 1631 400], 'InvertHardcopy', 'off');

% errorFcn = @(x) nanstd(x)/sqrt(size(x,1));
errorFcn = @(x) nanstd(x);
dims = {'x', 'y', 'z'};

difDifs = (modDifs - controlDifs);
mouseDifs = nan(length(mice),3,size(difDifs,3), 2); % [mouse X dim X kinDiffs X isBigStep]

for i = 1:length(mice)
    bins = strcmp({flat.mouse}, mice{i});
    mouseDifs(i,:,:,1) = nanmean(difDifs(bins & [flat.isBigStep],:,:),1);
    mouseDifs(i,:,:,2) = nanmean(difDifs(bins & ~[flat.isBigStep],:,:),1);
end

for i = 1:length(dims)
    
    % get difs for one dimension
    subplot(1,length(dims),i)    
    shadedErrorBar(times*1000, squeeze(mouseDifs(:,i,:,1)*1000), {@nanmean, errorFcn}, ...
        'lineprops', {'linewidth', 3, 'color', colors(1,:)}); hold on
    shadedErrorBar(times*1000, squeeze(mouseDifs(:,i,:,2)*1000), {@nanmean, errorFcn}, ...
        'lineprops', {'linewidth', 3, 'color', colors(2,:)});

    % pimp fig
    set(gca, 'xlim', xlims, 'ylim', ylims)
    ylabel(['\Delta' dims{i} ' (mm)'])
    xlabel('time (ms)')
    line([0 0], get(gca,'ylim'), 'color', [0 0 0])
end

blackenFig
% saveas(gcf, [getenv('OBSDATADIR') 'figures\reactionTimes\reactionTimes.png']);
% savefig(fullfile(getenv('OBSDATADIR'), 'figures\reactionTimes\reactionTimes.fig'))




%% plot single trial dif from control population

rows = 4;
cols = 4;
percentiles = [5 95];
trialInds = randperm(length(data), (rows*cols)/2);
xlims = [-.1 .2];
dim = 1;
    
figure('color', 'white');
errorFcn = @(x) 2*nanstd(x);

plotInd = 1;
for i = trialInds
    
    subaxis(rows, cols, plotInd, 'spacing', .02, 'margin', .08, 'padding', 0)
    if flat(i).isBigStep; colorInd=2; else; colorInd=1; end
    plot(times, squeeze(modLocations(i,dim,:)), 'color', colors(colorInd,:), 'linewidth', 3); hold on;
    shadedErrorBar(times, squeeze(controlLocations(i,dim,:,:))', {@nanmean, errorFcn}, ...
        'lineprops', {'linewidth', 3, 'color', [.5 .5 .5]});
    set(gca, 'box', 'off', 'xtick', times(1):.02:times(end), 'ytick', [], 'ycolor', 'white', 'xlim', xlims);
    plotInd = plotInd + 1;
    
    subaxis(rows, cols, plotInd, 'spacing', .02, 'margin', .08, 'padding', 0)
    plot(times, squeeze(controlLocations(i,dim,:,1)), 'color', [0 0 0], 'linewidth', 3); hold on;
    shadedErrorBar(times, squeeze(control2Locations(i,dim,:,:))', {@nanmean, errorFcn}, ...
        'lineprops', {'linewidth', 3, 'color', [.5 .5 .5]});
    set(gca, 'box', 'off', 'xtick', times(1):.02:times(end), 'ytick', [], 'ycolor', 'white', 'xlim', xlims);
    plotInd = plotInd + 1;
    
    
    line([0 0], get(gca,'ylim'), 'color', [0 0 0])
    
end
% blackenFig
pimpFig

% saveas(gcf, [getenv('OBSDATADIR') 'figures\reactionTimes\reactionTimesTrials.png']);



%% plot for each mouse separately



errorFcn = @(x) nanstd(x)/sqrt(size(x,1));
% errorFcn = @(x) nanstd(x);
dims = {'x', 'y', 'z'};

for mouse = 1:length(mice)
    figure('name', mice{mouse}, 'color', 'white', 'menubar', 'none', 'position', [2000 100*(mouse-1) 1000 200], 'InvertHardcopy', 'off');
    for i = 1:length(dims)
    
        % get diifs for one dimension
        bins = strcmp({flat.mouse}, mice{mouse});
        controlDifsSub = squeeze(controlDifs(bins,i,:));
        modDifsSub = squeeze(modDifs(bins,i,:));


        subplot(1,length(dims),i)

        shadedErrorBar(times*1000, controlDifsSub*1000, {@nanmean, errorFcn}, ...
            'lineprops', {'linewidth', 3, 'color', [0 0 0]}); hold on;
    %     shadedErrorBar(times*1000, -controlDifsSub*1000, {@nanmean, errorFcn}, ...
    %         'lineprops', {'linewidth', 3, 'color', [0 0 0]}); hold on;

        shadedErrorBar(times*1000, modDifsSub(~[flat(bins).isBigStep],:)*1000, {@nanmean, errorFcn}, ...
            'lineprops', {'linewidth', 3, 'color', colors(1,:)});
        shadedErrorBar(times*1000, modDifsSub([flat(bins).isBigStep],:)*1000, {@nanmean, errorFcn}, ...
            'lineprops', {'linewidth', 3, 'color', colors(2,:)});

        % pimp fig
        set(gca, 'xlim', xlims, 'ylim', [0 20])
        ylabel(['\Delta' dims{i} ' (mm)'])
        xlabel('time (ms)')
        line([0 0], get(gca,'ylim'), 'color', [0 0 0])
    end
    blackenFig
end




%% find and plot latency distributions
% !!! i think i computed the percentiles incorrectly // must check this!

% latencies = nan(1,length(data));
% for i = 1:length(data)
%     
%     modTrial = squeeze(modLocations(i,1,:));
%     controlErrs = prctile(squeeze(controlLocations(i,1,:,:))', percentiles, 1);
%     firstDeviatedInd = find(modTrial<controlErrs(1,:)' | modTrial>controlErrs(2,:)', 1, 'first');
%     if ~isempty(firstDeviatedInd)
%         latencies(i) = times(firstDeviatedInd);
%     else
%         latencies(i) = nan;
%     end
% end
% 
% figure; histogram(latencies, 20); pimpFig; blackenFig
% saveas(gcf, [getenv('OBSDATADIR') 'figures\reactionTimes\reactionTimeDistributions.png']);











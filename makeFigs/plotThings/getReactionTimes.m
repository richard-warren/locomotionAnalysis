% function getReactionTimes(data)

flat = getNestedStructFields(data, ...
    {'mouse', 'session', 'trial', 'preModPawKin', 'modPawKin', 'isBigStep', 'isLightOn', 'modPawContactInd', 'velAtWiskContact'});


%% settings
neighborNum = 40; % how many similar trials to use to construct control distribution
dt = .004; % 1/fps of camera
xlims = [-10 80]; % (ms)
ylims = [-20 20];
colors = [.25 1 1; .25 1 .25];
dimToMatch = 1;
minVel = .3;


% initializations
flat = flat([flat.modPawContactInd]>5 & [flat.modPawContactInd]<15); % add conditionals here
flat = flat(~[flat.isLightOn] & [flat.velAtWiskContact]>minVel); % add extra conditionals here // light on, velocity, only trials where there is a change in length?

swingMaxSmps = size(flat(1).modPawKin,2); % inherits max samples for swing locations from getKinematicData
times = linspace(-swingMaxSmps*dt, swingMaxSmps*dt, swingMaxSmps*2+1); % everything is centered around the moment of contact
contactInds = [flat.modPawContactInd];


% controlLocationsRaw1 = permute(cat(3, flat.preModPawKin), [3,1,2]);
controlLocationsRaw2 = permute(cat(3, flat.preModPawKin), [3,1,2]);
modLocationsRaw = permute(cat(3, flat.modPawKin), [3,1,2]);

% subtract x and y values s.t. they start at zero at the start of every trial
% controlLocationsRaw1 = controlLocationsRaw1 - repmat(controlLocationsRaw1(:,:,1), 1, 1, swingMaxSmps);
controlLocationsRaw2 = controlLocationsRaw2 - repmat(controlLocationsRaw2(:,:,1), 1, 1, swingMaxSmps);
modLocationsRaw = modLocationsRaw - repmat(modLocationsRaw(:,:,1), 1, 1, swingMaxSmps);


% initializations data containers
modLocations = nan(length(flat), 3, length(times));
controlLocations = nan(length(flat), 3, length(times), neighborNum);
modDifs = nan(length(flat), 3, length(times));
controlDifs = nan(length(flat), 3, length(times), neighborNum);


for i = 1:length(flat)
%     fprintf('\b\b\b\b\b')
%     fprintf('%.3f', i/length(flat))
    disp(i)
    
    % get temporal inds
    trialInds = (swingMaxSmps+1 : swingMaxSmps*2) - (contactInds(i)+1); % everything is centered at the momemnt of contact, s.t. moment of contact occurs at time 0

    % get inds for nearest neighbor trials
%     controlLocationsSub = permute(controlLocationsRaw2(:,:,1:contactInds(i)), [1,3,2]); % restrict to locations before moment of contact
%     controlLocationsSub = reshape(controlLocationsSub, length(flat), contactInds(i)*3); % concat xyz into a single dimension
%     modLocationsTrial = reshape(squeeze(modLocationsRaw(i,:,1:contactInds(i))), 1, contactInds(i)*3); % concat xyz into a single dimension
    
    controlLocationsSub = squeeze(controlLocationsRaw2(:,dimToMatch,1:contactInds(i)));
    modLocationsSub = squeeze(modLocationsRaw(i,dimToMatch,1:contactInds(i)))';
    
    inds2 = knnsearch(controlLocationsSub, modLocationsSub, 'k', neighborNum); % inds of nearest trials in terms of velocity
    
    controlLocations(i,:,trialInds,:) = permute(controlLocationsRaw2(inds2,:,:), [2,3,1]);

%     % get control locations and dif of control locations from avg control locations
%     for j = 1:length(inds2)
% 
%         % get single control locations
%         controlLocations(i,:,trialInds,j) = controlLocationsRaw2(inds2(j),:,:);
% 
%         inds1 = knnsearch(controlVels1', controlVels2(inds2(j)), 'k', neighborNum);
% 
%         % get dif between average of remaining control locations and control trial
%         controlMean = nanmean(controlLocationsRaw1(inds1,:,:), 1);
%         controlDifs(i,:,trialInds,j) = controlLocations(i,:,trialInds,j) - controlMean;
%     end

    % get mod locations and dif from control
    modLocations(i,:,trialInds) = modLocationsRaw(i,:,:);

    % get dif between average of control locations and mod locations
    controlMean = nanmean(controlLocationsRaw2(inds2,:,:), 1);
    modDifs(i,:,trialInds) = modLocations(i,:,trialInds) - controlMean;
end




figure('color', 'white', 'menubar', 'none', 'position', [2052 244 1631 491], 'InvertHardcopy', 'off');

% errorFcn = @(x) nanstd(x)/sqrt(size(x,1));
errorFcn = @(x) nanstd(x);
dims = {'x', 'y', 'z'};

for i = 1:length(dims)
    
    % get diifs for one dimension
%     controlDifsSub = squeeze(controlDifs(:,i,:,:));
%     controlDifsSub = reshape(permute(controlDifsSub, [1 3 2]), [], length(times));
    modDifsSub = squeeze(modDifs(:,i,:));
    
    
    subplot(1,length(dims),i)
    
%     shadedErrorBar(times*1000, controlDifsSub*1000, {@nanmean, errorFcn}, ...
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






%% plot single trial dif from control population

rows = 4;
cols = 3;
percentiles = [5 95];
trialInds = randperm(length(data), (rows*cols));
xlims = [-.1 .2];
dim = 3;
    
figure('color', 'white');

for i = 1:(rows*cols)
    subaxis(rows, cols, i, 'spacing', .02, 'margin', .08, 'padding', 0)
    
    if flat(trialInds(i)).isBigStep; colorInd=2; else; colorInd=1; end
    plot(times, squeeze(modLocations(trialInds(i),dim,:)), 'color', colors(colorInd,:), 'linewidth', 3); hold on;
    
%     controlMean = nanmean(squeeze(controlLocations(trialInds(i),dim,:,:)), 2)';
%     controlErrs = prctile(squeeze(controlLocations(trialInds(i),dim,:,:))', percentiles, 1);
%     controlErrs = abs(controlErrs - repmat(controlMean,2,1));
%     shadedErrorBar(times, controlMean, controlErrs, 'lineprops', {'linewidth', 3, 'color', [0 0 0]});
    shadedErrorBar(times, squeeze(controlLocations(trialInds(i),dim,:,:))', {@nanmean, errorFcn}, ...
        'lineprops', {'linewidth', 3, 'color', [0 0 0]});
    
    set(gca, 'box', 'off', 'xtick', times(1):.02:times(end), 'ytick', [], 'ycolor', 'white', 'xlim', xlims);
    line([0 0], get(gca,'ylim'), 'color', [0 0 0])
    
end
% blackenFig
pimpFig

% saveas(gcf, [getenv('OBSDATADIR') 'figures\reactionTimes\reactionTimesTrials.png']);



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











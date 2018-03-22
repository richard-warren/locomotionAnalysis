

% settings
sessions = {'180122_001', '180122_002', '180122_003', ...
            '180123_001', '180123_002', '180123_003', ...
            '180124_001', '180124_002', '180124_003', ...
            '180125_001', '180125_002', '180125_003'};
neighborNum = 40; % how many similar trials to use to construct control distribution
dt = .004; % 1/fps of camera
dtInterp = .001;


% initializations
dataRaw = getKinematicData(sessions);
data = dataRaw([dataRaw.oneSwingOneStance]);
swingMaxSmps = size(data(1).modifiedLocations{1},3); % inherits max samples for swing locations from getKinematicData
times = linspace(-swingMaxSmps*dt, swingMaxSmps*dt, swingMaxSmps*2+1);
timesInterp = times(1):dtInterp:times(end);

%%


controlVels1 = cellfun(@(x) x(1,3), {data.controlWheelVels});
controlVels2 = cellfun(@(x) x(2,3), {data.controlWheelVels});
modVels = cellfun(@(x) x(1,3), {data.modifiedWheelVels});
contactInds = reshape([data.pawObsPosInd]',4,length(data))'; contactInds = contactInds(:,3);
contactInds = contactInds + 1; % !!! this is a hack - need to fix getKinematicData s.t. no inds are zero

controlLocationsRaw1 = cellfun(@(x) x{3}(end,:,:), {data.controlLocations}, 'uniformoutput', 0);
controlLocationsRaw1 = cat(1,controlLocationsRaw1{:});
controlLocationsRaw2 = cellfun(@(x) x{3}(end,:,:), {data.controlLocations}, 'uniformoutput', 0);
controlLocationsRaw2 = cat(1,controlLocationsRaw2{:});
modLocationsRaw = cellfun(@(x) x{3}(1,:,:), {data.modifiedLocations}, 'uniformoutput', 0);
modLocationsRaw = cat(1,modLocationsRaw{:});

% !!! for each swing, fill in final NaNs with last x,y values to account for stance position
% NOTE: this is a hack that will not work if eventually you want to use x y and z for euclidian distance
% if would be better to store stride locations AND subsequent locations in getKinematicData, perhaps along with a thing that tells you stride duration, so you know which inds correspond to the actual stride
% then you don't need to infer what the values WOULD have been // instead you have them directly
for i = 1:length(data)
    
    nanInd = find(isnan(controlLocationsRaw1(i,1,:)),1,'first');
    controlLocationsRaw1(i,1,nanInd:end) = controlLocationsRaw1(i,1,nanInd-1);
    
    nanInd = find(isnan(controlLocationsRaw2(i,1,:)),1,'first');
    controlLocationsRaw2(i,1,nanInd:end) = controlLocationsRaw2(i,1,nanInd-1);
    
    nanInd = find(isnan(modLocationsRaw(i,1,:)),1,'first');
    modLocationsRaw(i,1,nanInd:end) = modLocationsRaw(i,1,nanInd-1);
    
end

% subtract x and y values s.t. they start at zero at the start of every trial
controlLocationsRaw1 = controlLocationsRaw1 - repmat(controlLocationsRaw1(:,:,1), 1, 1, swingMaxSmps);
controlLocationsRaw2 = controlLocationsRaw2 - repmat(controlLocationsRaw2(:,:,1), 1, 1, swingMaxSmps);
modLocationsRaw = modLocationsRaw - repmat(modLocationsRaw(:,:,1), 1, 1, swingMaxSmps);



% initializations data containers
modLocations = nan(length(data), 2, length(times));
controlLocations = nan(length(data), 2, length(times), neighborNum);
modDifs = nan(length(data), 2, length(times));
controlDifs = nan(length(data), 2, length(times), neighborNum);


for i = 1:length(data)
    
    % get temporal inds
    trialInds = (swingMaxSmps+1 : swingMaxSmps*2) - (contactInds(i)+1);
    
    % get inds for nearest neighbor trials
    inds2 = knnsearch(controlVels2', modVels(i), 'k', neighborNum);
    
    % get control locations and dif of control locations from avg control locations
    for j = 1:length(inds2)
        
        % get single control locations
        controlLocations(i,:,trialInds,j) = controlLocationsRaw2(inds2(j),:,:);
        
        inds1 = knnsearch(controlVels1', controlVels2(inds2(j)), 'k', neighborNum);
        
        % get dif between average of remaining control locations and control trial
        controlMean = nanmean(controlLocationsRaw1(inds1,:,:), 1);
        controlDifs(i,:,trialInds,j) = (controlLocations(i,:,trialInds,j) - controlMean);
    end
    
    % get mod locations and dif from control
    modLocations(i,:,trialInds) = modLocationsRaw(i,:,:);
    
    % get dif between average of control locations and mod locations
    controlMean = nanmean(controlLocationsRaw2(inds2,:,:), 1);
    modDifs(i,:,trialInds) = (modLocations(i,:,trialInds) - controlMean);    
end


% get x values
controlDifsX = squeeze(controlDifs(:,1,:,:));
controlDifsX = reshape(permute(controlDifsX, [1 3 2]), [], length(times));
modDifsX = squeeze(modDifs(:,1,:));

% interpolate
% controlDifsXInterp = interpWithNans(controlDifsX, times, timesInterp, 'pchip');
% modDifsXInterp = interpWithNans(modDifsX, times, timesInterp, 'pchip');


% plot mean difs
xlims = [-10 100];
ylims = [-20 30];
colors = winter(2);

close all; figure('color', 'white', 'menubar', 'none', 'position', [680   215   528   763]);
numModSteps = reshape([data.modStepNum],4,length(data))';
fcns = {@(x) nanstd(x)/sqrt(size(x,1)), @(x) nanstd(x)};

for i = 1:length(fcns)
    
    subplot(2,1,i)
    
    shadedErrorBar(times*1000, controlDifsX*1000, {@nanmean, fcns{i}}, ...
        'lineprops', {'linewidth', 3, 'color', [0 0 0]}); hold on;
    shadedErrorBar(times*1000, modDifsX(numModSteps(:,3)>1,:)*1000, {@nanmean, fcns{i}}, ...
        'lineprops', {'linewidth', 3, 'color', colors(1,:)});
    shadedErrorBar(times*1000, modDifsX(numModSteps(:,3)==1,:)*1000, {@nanmean, fcns{i}}, ...
        'lineprops', {'linewidth', 3, 'color', colors(2,:)});

    % pimp fig
    set(gca, 'xlim', xlims, 'ylim', ylims)
    ylabel('\Deltax (mm)')
    xlabel('time (ms)')
    line([0 0], get(gca,'ylim'), 'color', [0 0 0])
end

saveas(gcf, [getenv('OBSDATADIR') 'figures\reactionTimes.png']);


%% plot single trial dif from control population

rows = 4;
cols = 3;
percentiles = [5 95];
trialInds = randperm(length(data), (rows*cols));
xlims = [-.1 .2];

close all; figure('color', 'white');

for i = 1:(rows*cols)
    subaxis(rows, cols, i, 'spacing', .02, 'margin', .08, 'padding', 0)
    
    if numModSteps(trialInds(i),3)==1; colorInd=2; else; colorInd=1; end
    plot(times, squeeze(modLocations(trialInds(i),1,:)), 'color', colors(colorInd,:), 'linewidth', 3); hold on;
    
    controlMean = nanmean(squeeze(controlLocations(trialInds(i),1,:,:)), 2)';
    controlErrs = prctile(squeeze(controlLocations(trialInds(i),1,:,:))', percentiles, 1);
    controlErrs = abs(controlErrs - repmat(controlMean,2,1));
    shadedErrorBar(times, controlMean, controlErrs, 'lineprops', {'linewidth', 3, 'color', [0 0 0]});
    
    set(gca, 'box', 'off', 'ytick', [], 'xtick', [], 'xcolor', 'white', 'ycolor', 'white', 'xlim', xlims);
    line([0 0], get(gca,'ylim'), 'color', [0 0 0])
    
end

pimpFig

saveas(gcf, [getenv('OBSDATADIR') 'figures\reactionTimesTrials.png']);


%% find and plot latency distributions

latencies = 1:length(data);

for i = 1:length(data)
    
    modTrial = squeeze(modLocations(i,1,:));
    controlErrs = prctile(squeeze(controlLocations(i,1,:,:))', percentiles, 1);
    firstDeviatedInd = find(modTrial<controlErrs(1,:)' | modTrial>controlErrs(2,:)', 1, 'first');
    if ~isempty(firstDeviatedInd)
        latencies(i) = times(firstDeviatedInd);
    else
        latencies(i) = nan;
    end
end

figure; histogram(latencies(latencies>0));



























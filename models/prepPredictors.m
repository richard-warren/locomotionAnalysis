function prepPredictors(session)

% settings
s.dt = .002;
s.velTime = .05;  % (s) velocity is computed over this interval
s.percentileLims = [.1 99.9];  % remove and interpolate tracking outside this percentile range
s.plotPredictors = true;

% initializations
% if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end % reassign settings passed in varargin
fprintf('%s: preparing predictors...\n', session)
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'frameTimeStamps', 'frameTimeStampsWisk', 'wheelPositions', 'wheelTimes', ...
    'bodyAngles', 'whiskerAngle', 'pixelsPerM', 'wheelCenter', 'wheelRadius', ...
    'lickTimes', 'touchesPerPaw', 'touchClassNames', 'obsOnTimes', 'obsOffTimes', ...
    'rewardTimes', 'isRewardSurprise', 'omissionTimes', 'wiskContactTimes', ...
    'obsLightOnTimes', 'obsLightOffTimes');
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'run.mat'), 'obsLight');
locationsRun = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw.csv'));
locationsWisk = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw_wisk.csv'));
pawNames = {'paw1LH', 'paw2LF', 'paw3RF', 'paw4RH'};
dimensions = {'x', 'y', 'z'};
fps = 1/nanmedian(diff(frameTimeStamps));
scoreThresh = getScoreThresh(session, 'trackedFeaturesRaw_metadata.mat');  % scoreThresh depends on whether deeplabcut (old version) or deepposekit was used
validTimes = ~isnan(frameTimeStamps);
validTimesWisk = ~isnan(frameTimeStampsWisk);




% ----------
% CONTINUOUS
% ----------

% define time grid and empty table
tMin = max([frameTimeStamps(1), frameTimeStampsWisk(1), 0]);  % start with recording device that starts last (spike starts at zero)
tMax = min([frameTimeStamps(end), frameTimeStampsWisk(end), wheelTimes(end)]);  % start with recording device that ends first
t = tMin:s.dt:tMax;
cntPredictors = table(t', 'VariableNames', {'t'});

% wheel velocity
vel = getVelocity(wheelPositions, s.velTime, 1/nanmedian(diff(wheelTimes)));
vel = interp1(wheelTimes, vel, t);
cntPredictors.velocity = vel(:);

% paw locations and velocity
[pawXYZ, pawXYZ_pixels] = getPawXYZ(session);
for i = 1:length(pawNames)
    for j = 1:3
        pos = pawXYZ.(pawNames{i})(:,j);
        pos = interp1(frameTimeStamps(validTimes), pos(validTimes), t);
        cntPredictors.([pawNames{i} '_' dimensions{j}]) = pos(:);
    end
end

% body angle
cntPredictors.bodyAngle = interp1(frameTimeStamps(validTimes), bodyAngles(validTimes), t)';

% whisker angle
cntPredictors.whiskerAngle = interp1(frameTimeStampsWisk(validTimesWisk), whiskerAngle(validTimesWisk), t)';

% butt height
[buttHeight, buttConfidence] = deal(locationsRun.tailBase_top_1, locationsRun.tailBase_top_2);
buttHeight(buttConfidence<scoreThresh) = nan;
buttHeight = fillmissing(buttHeight, 'pchip');
buttHeight = ((wheelCenter(2)-wheelRadius) - buttHeight) / pixelsPerM; % flip z and set s.t. top of wheel is zero
buttHeight = interp1(frameTimeStamps(validTimes), buttHeight(validTimes), t);
cntPredictors.buttHeight = buttHeight(:);

% jaw (first pc projection)
[jawX, jawZ, jawConfidence] = deal(locationsWisk.jaw, locationsWisk.jaw_1, locationsWisk.jaw_2);
jaw = pcProject([jawX, jawZ], jawConfidence, true);
jaw = interp1(frameTimeStampsWisk(validTimesWisk), jaw(validTimesWisk), t);
cntPredictors.jaw = jaw(:);

% ear (first pc projection)
[earX, earZ, earConfidence] = deal(locationsRun.ear, locationsRun.ear_1, locationsRun.ear_2);
ear = pcProject([earX, earZ], earConfidence, true);
ear = interp1(frameTimeStamps(validTimes), ear(validTimes), t);
cntPredictors.ear = ear(:);

% nose (first pc projection)
[noseX, noseZ, noseConfidence] = deal(locationsWisk.nose, locationsWisk.nose_1, locationsWisk.nose_2);
nose = pcProject([noseX, noseZ], noseConfidence, true);
nose = interp1(frameTimeStampsWisk(validTimesWisk), nose(validTimesWisk), t);
cntPredictors.nose = nose(:);

% % whisker pad (first pc projection)
% [padX, padZ, padConfidence] = deal(locationsWisk.wisk_pad, locationsWisk.wisk_pad_1, locationsWisk.wisk_pad_2);
% pad = pcProject([padX, padZ], padConfidence, true);
% pad = interp1(frameTimeStampsWisk(validTimesWisk), pad(validTimesWisk), t);
% cntPredictors.pad = pad(:);

% satiation
lickBins = histcounts(lickTimes, [t (t(end)+s.dt)]);
satiation = cumsum(lickBins) / sum(lickBins);
cntPredictors.satiation = satiation(:);

% velocity for continuous predictors
names = cntPredictors.Properties.VariableNames;
exclude = {'t', 'velocity', 'satiation'};  % don't compute velocity for these predictors
for i = 1:length(names)
    if ~ismember(names{i}, exclude)
        vel = getVelocity(cntPredictors.(names{i}), s.velTime, 1/s.dt);
        cntPredictors.([names{i} '_vel']) = vel(:);
    end
end




% -----
% EPOCH
% -----

% swing/stance
xzLocations = nan(length(frameTimeStamps), 2, 4);  % (time X xz X paw)
for i = 1:length(pawNames); xzLocations(:,:,i) = pawXYZ_pixels.(pawNames{i})(:,[1,3]); end  % format for getStanceBins
stanceBins = getStanceBins(frameTimeStamps, xzLocations, wheelPositions, wheelTimes, wheelCenter, wheelRadius, fps, pixelsPerM);
for i = 1:4
    epochPredictors.([pawNames{i} '_stance']) = logicalToEpochs(stanceBins(:,i), frameTimeStamps);
end

% obstacle
epochPredictors.obstacle = [obsOnTimes(:) obsOffTimes(:)];

% obstacle light
epochPredictors.light = [obsLightOnTimes(:) obsLightOffTimes(:)];  % assumes same number of on and off times

% reward
epochPredictors.reward = [rewardTimes(1:end-1), rewardTimes(2:end)];

% stride (start of swing to start of next swing)
for i = 1:4
    swingStartTimes = frameTimeStamps(find(diff(stanceBins(:,i))==-1)+1);
    epochPredictors.([pawNames{i} '_stride']) = [swingStartTimes(1:end-1) swingStartTimes(2:end)];
end




% ------
% EVENTS
% ------

% whisker contact
eventPredictors.whiskerContact = wiskContactTimes;

% licks
eventPredictors.lick = lickTimes;

% rewards
eventPredictors.reward_all = sort([rewardTimes; omissionTimes]);
eventPredictors.reward_normal = rewardTimes(~isRewardSurprise);
eventPredictors.reward_omission = omissionTimes;
eventPredictors.reward_surprise = rewardTimes(isRewardSurprise);

% paw contacts
medianFiltering = 5;  % (samples) window size for median filtering to remove 'blips'
minDif = 1;  % (seconds) paw contacts must be at least this far apart in time
touchesPerPaw(touchesPerPaw==0) = find(strcmp(touchClassNames, 'no_touch'));  % recode no_touch class
include = {'fore_dorsal', 'fore_ventral', 'hind_dorsal', 'hind_ventral_low'};
for i = 1:4
    for j = {'dorsal', 'ventral'}
        touchTypes = touchClassNames(touchesPerPaw(:,i));
        touchBins = contains(touchTypes, j{1}) & ismember(touchTypes, include);
        touchBins = logical(medfilt1(single(touchBins), medianFiltering));
        touchTimes = frameTimeStamps(find(diff(touchBins)==1)+1);
        if ~isempty(touchTimes)
            touchTimes = touchTimes([true; diff(touchTimes)>minDif]);  % remove touches that are too close together in time
        end
        eventPredictors.([pawNames{i} '_contact_' j{1}]) = touchTimes;
    end
end




% ------------------
% PLOT YOUR FACE OFF
% ------------------

% make nice little plot
if s.plotPredictors
    
    % settings
    dz = 6;  % (standard deviation) vertical separation
    xWidth = 40;
    xLims = [0 xWidth] + randi(round(tMax-xWidth*2));
    ommit = {'paw1LH_stride', 'paw2LF_stride', 'paw3RF_stride', 'paw4RH_stride', ...
        'reward', 'reward_surprise', 'reward_omission', 'reward_normal'};
    
    % initialiations
    set(0, 'DefaultAxesTickLabelInterpreter', 'none')
    figure('color', 'white', 'position', [125.00 93.00 1666.00 886.00]); hold on
    eventNames = fieldnames(eventPredictors); eventNames = eventNames(~ismember(eventNames, ommit));
    epochNames = fieldnames(epochPredictors); epochNames = epochNames(~ismember(epochNames, ommit));
    cntNames = cntPredictors.Properties.VariableNames(2:end); cntNames = cntNames(~ismember(cntNames, ommit));
    allNames = [epochNames; cntNames'; eventNames];
    rows = length(allNames);
    colors = lines(rows);
    y = 0:dz:dz*rows-1;
    
    % epoch
    plot(xLims, repmat(y(1:length(epochNames)),2,1)', 'color', [.4 .4 .4])  % horizontal lines for each epoch predictor
    for i = 1:length(epochNames)
        epoch = epochPredictors.(epochNames{i});
        epoch = epoch(any(epoch>xLims(1) & epoch<xLims(2),2),:);  % only plot epochs within xLims
        plot(epoch', repelem(y(i), 2), 'LineWidth', 5, 'color', colors(i,:))
    end
    
    % continuous
    dy = length(epochNames);
    bins = t>xLims(1) & t<xLims(2);
    yOffsets = y((1:length(cntNames))+dy);
    plot(t(bins), zscore(table2array(cntPredictors(bins,2:end))) + yOffsets);
    
    % events
    dy = length(epochNames) + length(cntNames);
    plot(xLims, repmat(y((1:length(eventNames))+dy),2,1)', 'color', [.4 .4 .4])  % horizontal lines for each event predictor
    for i = 1:length(eventNames)
        x = eventPredictors.(eventNames{i});
        x = x(x>xLims(1) & x<xLims(2));
        scatter(x, repelem(y(i+dy), length(x)), 20, colors(i,:), 'filled')
    end
    
    % fancify
    set(gca, 'xlim', xLims, 'ytick', y, 'YTickLabel', allNames, 'YLim', [y(1) y(end)], 'TickDir', 'out')
end



% ----
% SAVE
% ----

dirName = fullfile(getenv('OBSDATADIR'), 'sessions', session, 'modelling');
if ~exist(dirName, 'dir'); mkdir(dirName); end
save(fullfile(dirName, 'predictors.mat'), 'cntPredictors', 'epochPredictors', 'eventPredictors')
if s.plotPredictors; savefig(gcf, fullfile(dirName, 'predictors')); end




% ---------
% FUNCTIONS
% ---------

function projection = pcProject(sig, confidence, normalize)
    % project signal onto first principal component
    % signal is (time X dimension) matrix
    
    if ~exist('normalize', 'var'); normalize = false; end
    
    sig(confidence<scoreThresh,:) = nan;
    pcs = pca(sig);
    projection = sig * pcs(:,1);  % project onto first pc
    lims = prctile(projection, s.percentileLims);
    projection(projection<lims(1) | projection>lims(2)) = nan;
    projection = fillmissing(projection, 'pchip');
    
    if normalize
        projection = (projection-nanmean(projection)) / pixelsPerM;
    end
end


function epochs = logicalToEpochs(logical, logicalTimes)
    % give binary signal, produces (n X 2) matrix where each row is the
    % start and end time for epochs where the signal is true
    startInds = find(diff(logical)==1)+1;
    endInds = find(diff(logical)==-1)+1;
    
    % pad to make sure starts with start and ends with end
    if startInds(1)>endInds(1); startInds = [1; startInds]; end
    if endInds(end)<startInds(end); endInds = [endInds; length(logicalTimes)]; end
    
    epochs = [logicalTimes(startInds), logicalTimes(endInds)];
    epochs = epochs(~any(isnan(epochs),2),:);  % remove rows with nans
end

end




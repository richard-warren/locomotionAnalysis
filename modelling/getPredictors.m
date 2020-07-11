function getPredictors(session, varargin)

% prepares table containing all predictors used for neural encoding models
% // data are either events, epochs (e.g. isObsOn), or continuous // if
% predictor is unusable, 'include' is set to 0


% settings
s.dt = .002;  % (s) everything interpolated onto new time axis with this temporal resolution
s.velTime = .05;  % (s) velocity is computed over this interval
s.percentileLims = [.1 99.9];  % remove and interpolate tracking outside this percentile range
s.plotPredictors = true;



% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end % reassign settings passed in varargin
fprintf('%s: preparing predictors...\n', session)
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'frameTimeStamps', 'frameTimeStampsWisk', 'wheelPositions', 'wheelTimes', ...
    'bodyAngles', 'whiskerAngle', 'pixelsPerM', 'wheelCenter', 'wheelRadius', ...
    'lickTimes', 'touchesPerPaw', 'touchClassNames', 'obsOnTimes', 'obsOffTimes', ...
    'rewardTimes', 'isRewardSurprise', 'omissionTimes', 'wiskContactTimes', ...
    'obsLightOnTimes', 'obsLightOffTimes');
locationsRun = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw.csv'));
locationsWisk = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw_wisk.csv'));
pawNames = {'paw1LH', 'paw2LF', 'paw3RF', 'paw4RH'};
dimensions = {'x', 'y', 'z'};
fps = 1/nanmedian(diff(frameTimeStamps));
scoreThresh = getScoreThresh(session, 'trackedFeaturesRaw_metadata.mat');  % scoreThresh depends on whether deeplabcut (old version) or deepposekit was used
validTimes = ~isnan(frameTimeStamps);
validTimesWisk = ~isnan(frameTimeStampsWisk);
predictors = table({}, categorical({}, {'event', 'epoch', 'continuous'}), {}, [], ...
    'VariableNames', {'data', 'type', 't', 'include'}, 'RowNames', {});



% ----------
% CONTINUOUS
% ----------

% define time grid and empty table
tMin = max([frameTimeStamps(1), frameTimeStampsWisk(1), 0]);  % start with recording device that starts last (spike starts at zero)
tMax = min([frameTimeStamps(end), frameTimeStampsWisk(end), wheelTimes(end)]);  % start with recording device that ends first
t = tMin:s.dt:tMax;

% wheel velocity
vel = getVelocity(wheelPositions, s.velTime, 1/nanmedian(diff(wheelTimes)));
vel = interp1(wheelTimes, vel, t);
addPredictor('velocity', vel(:)', 'continuous', t)

% paw position and pahse
[pawXYZ, pawXYZ_pixels] = getPawXYZ(session);
for i = 1:length(pawNames)
    for j = 1:3
        % position
        pos = pawXYZ.(pawNames{i})(:,j);
        pos = interp1(frameTimeStamps(validTimes), pos(validTimes), t);
        addPredictor([pawNames{i} '_' dimensions{j}], pos, 'continuous', t)
        
        % phase
        if j==1
            pos = highpass(pos, .1, 1/median(diff(t)));
            phase = angle(hilbert(pos));
            addPredictor([pawNames{i} '_phase'], phase, 'continuous', t)
        end
    end
end

% body angle
bodyAngle = interp1(frameTimeStamps(validTimes), bodyAngles(validTimes), t)';
addPredictor('bodyAngle', bodyAngle, 'continuous', t)

% whisker angle
whiskerAngle = interp1(frameTimeStampsWisk(validTimesWisk), whiskerAngle(validTimesWisk), t)';
addPredictor('whiskerAngle', whiskerAngle, 'continuous', t)

% butt height
[buttHeight, buttConfidence] = deal(locationsRun.tailBase_top_1, locationsRun.tailBase_top_2);
buttHeight(buttConfidence<scoreThresh) = nan;
buttHeight = fillmissing(buttHeight, 'pchip');
buttHeight = ((wheelCenter(2)-wheelRadius) - buttHeight) / pixelsPerM; % flip z and set s.t. top of wheel is zero
buttHeight = interp1(frameTimeStamps(validTimes), buttHeight(validTimes), t);
addPredictor('buttHeight', buttHeight, 'continuous', t)

% jaw (first pc projection)
[jawX, jawZ, jawConfidence] = deal(locationsWisk.jaw, locationsWisk.jaw_1, locationsWisk.jaw_2);
jaw = pcProject([jawX, jawZ], jawConfidence, true, 'jaw');
jaw = interp1(frameTimeStampsWisk(validTimesWisk), jaw(validTimesWisk), t);
addPredictor('jaw', jaw, 'continuous', t)

% ear (first pc projection)
[earX, earZ, earConfidence] = deal(locationsRun.ear, locationsRun.ear_1, locationsRun.ear_2);
ear = pcProject([earX, earZ], earConfidence, true, 'ear');
ear = interp1(frameTimeStamps(validTimes), ear(validTimes), t);
addPredictor('ear', ear, 'continuous', t)

% nose (first pc projection)
[noseX, noseZ, noseConfidence] = deal(locationsWisk.nose, locationsWisk.nose_1, locationsWisk.nose_2);
nose = pcProject([noseX, noseZ], noseConfidence, true, 'nose');
nose = interp1(frameTimeStampsWisk(validTimesWisk), nose(validTimesWisk), t);
addPredictor('nose', nose, 'continuous', t)

% satiation
lickBins = histcounts(lickTimes, [t (t(end)+s.dt)]);
satiation = cumsum(lickBins) / sum(lickBins);
addPredictor('satiation', satiation, 'continuous', t)

% velocity for continuous predictors
% (expect for those listed in 'exclude')
excludeVars = {'velocity', 'satiation', 'paw1LH_phase', 'paw2LF_phase', 'paw3RF_phase', 'paw4RH_phase'};  % don't compute velocity for these predictors
for row = predictors.Properties.RowNames'
    if ~ismember(row{1}, excludeVars)
        vel = getVelocity(predictors{row{1}, 'data'}{1}, s.velTime, 1/s.dt);
        addPredictor([row{1} '_vel'], vel, 'continuous', t)
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
    epoch = logicalToEpochs(stanceBins(:,i), frameTimeStamps);
    addPredictor([pawNames{i} '_stance'], epoch, 'epoch')
end

% obstacle
addPredictor('obstacle', [obsOnTimes(:) obsOffTimes(:)], 'epoch')

% obstacle no light
addPredictor('light', [obsLightOnTimes(:) obsLightOffTimes(:)], 'epoch')

% reward
addPredictor('reward', [rewardTimes(1:end-1), rewardTimes(2:end)], 'epoch')

% stride (start of swing to start of next swing)
for i = 1:4
    swingStartTimes = frameTimeStamps(find(diff(stanceBins(:,i))==-1)+1);
    addPredictor([pawNames{i} '_stride'], [swingStartTimes(1:end-1) swingStartTimes(2:end)], 'epoch')
end



% ------
% EVENTS
% ------

% whisker contact
addPredictor('whiskerContact', wiskContactTimes, 'event')

% licks
addPredictor('lick', lickTimes, 'event')

% rewards
addPredictor('reward_all', sort([rewardTimes; omissionTimes]), 'event')
addPredictor('reward_normal', rewardTimes(~isRewardSurprise), 'event')
addPredictor('reward_omission', omissionTimes, 'event')
addPredictor('reward_surprise', rewardTimes(isRewardSurprise), 'event')

% paw contacts
medianFiltering = 5;  % (samples) window size for median filtering to remove 'blips'
minDif = 1;  % (seconds) paw contacts must be at least this far apart in time
touchesPerPaw(touchesPerPaw==0) = find(strcmp(touchClassNames, 'no_touch'));  % recode no_touch class
includeTouches = {'fore_dorsal', 'fore_ventral', 'hind_dorsal', 'hind_ventral_low'};
for i = 1:4
    for j = {'dorsal', 'ventral'}
        touchTypes = touchClassNames(touchesPerPaw(:,i));
        touchBins = contains(touchTypes, j{1}) & ismember(touchTypes, includeTouches);
        touchBins = logical(medfilt1(single(touchBins), medianFiltering));
        touchTimes = frameTimeStamps(find(diff(touchBins)==1)+1);
        if ~isempty(touchTimes)
            touchTimes = touchTimes([true; diff(touchTimes)>minDif]);  % remove touches that are too close together in time
        end
        addPredictor([pawNames{i} '_contact_' j{1}], touchTimes, 'event')
    end
end




% -----
% PLOT!
% -----

% make nice little plot
if s.plotPredictors
    
    % settings
    dz = 6;  % (standard deviation) vertical separation
    xWidth = 40;  % (seconds) range of x axis
    xLims = [0 xWidth] + randi(round(tMax-xWidth*2));
    ommit = {'paw1LH_stride', 'paw2LF_stride', 'paw3RF_stride', 'paw4RH_stride', ...
        'reward', 'reward_surprise', 'reward_omission', 'reward_normal'};
    
    % initialiations
    set(0, 'DefaultAxesTickLabelInterpreter', 'none')
    figure('name', session, 'color', 'white', 'position', [125.00 93.00 1666.00 886.00]); hold on
    allInds = find(~ismember(predictors.Properties.RowNames, ommit));
    colors = lines(length(allInds));
    y = 0:dz:dz*length(allInds)-1;
    
    % epoch
    epochInds = allInds(predictors.type(allInds)=='epoch');
    plot(xLims, repmat(y(1:length(epochInds)),2,1)', 'color', [.4 .4 .4])  % horizontal lines for each epoch predictor
    for i = 1:length(epochInds)
        epoch = predictors{epochInds(i),'data'}{1};
        epoch = epoch(any(epoch>xLims(1) & epoch<xLims(2),2),:);  % only plot epochs within xLims
        if ~isempty(epoch)
            plot(epoch', repelem(y(i), 2), 'LineWidth', 5, 'color', colors(i,:))
        end
    end
    
    % continuous
    contInds = allInds(predictors.type(allInds)=='continuous');
    bins = t>xLims(1) & t<xLims(2);
    yOffsets = y((1:length(contInds))+length(epochInds));
    allCont = cat(1,predictors{predictors.type=='continuous','data'}{:});
    plot(t(bins), zscore(allCont(:,bins)') + yOffsets);
    
    % events
    eventInds = allInds(predictors.type(allInds)=='event');
    dy = length(epochInds) + length(contInds);
    plot(xLims, repmat(y((1:length(eventInds))+dy),2,1)', 'color', [.4 .4 .4])  % horizontal lines for each event predictor
    for i = 1:length(eventInds)
        x = predictors{eventInds(i),'data'}{1};
        x = x(x>xLims(1) & x<xLims(2));
        scatter(x, repelem(y(i+dy), length(x)), 20, colors(i+dy,:), 'filled')
    end
    
    % fancify
    set(gca, 'xlim', xLims, 'ytick', y, ...
        'YTickLabel', predictors.Properties.RowNames([epochInds; contInds; eventInds]), ...
        'YLim', [y(1) y(end)], 'TickDir', 'out')
end



% ----
% SAVE
% ----

dirName = fullfile(getenv('OBSDATADIR'), 'sessions', session, 'modelling');
if ~exist(dirName, 'dir'); mkdir(dirName); end
save(fullfile(dirName, 'predictors.mat'), 'predictors')
if s.plotPredictors
    saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures', 'modelling', 'predictors', sprintf('%s predictors.png', session)));
end



% ---------
% FUNCTIONS
% ---------

function projection = pcProject(sig, confidence, normalize, name)
    % project signal onto first principal component
    % signal is (time X dimension) matrix
    
    validBins = confidence>=scoreThresh;
    if mean(validBins)>.5  % at least half of the frame must be tracked
        sig(~validBins,:) = nan;
        pcs = pca(sig);
        projection = sig * pcs(:,1);  % project onto first pc
        lims = prctile(projection, s.percentileLims);
        projection(projection<lims(1) | projection>lims(2)) = nan;
        projection = fillmissing(projection, 'pchip');
    
        if normalize
            projection = (projection-nanmean(projection)) / pixelsPerM;
        end
    else
        projection = nan(1,length(confidence));
        fprintf('  %s: %s tracking unsuccessful\n', session, name)
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


function addPredictor(name, data, type, t)
    % extends predictors table by add a new row
    % include=0 when there are no events, no epochs, or all of the
    % continuous dta are NaN (occurs when pcProject fails)
    
    if ~exist('include', 'var'); include = true; end
    if ~exist('t', 'var'); t = []; end
    
    % ensure homogeneous orientation
    if strcmp(type, 'event')
        data=data(:);
        include = length(data)>1;
    elseif strcmp(type, 'epoch')
        include = size(data,1)>1;
    elseif strcmp(type, 'continuous')
        data=data(:)';
        include = ~all(isnan(data));
    end
    
    newRow = table({data}, categorical({type}, {'event', 'epoch', 'continuous'}), {t}, include, ...
        'VariableNames', {'data', 'type', 't', 'include'}, 'RowNames', {name});
    predictors = [predictors; newRow];
end

end




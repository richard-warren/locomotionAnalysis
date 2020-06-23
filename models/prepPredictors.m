function prepPredictors(session)

% settings
s.dt = .002;
s.velTime = .05;  % (s) velocity is computed over this interval
s.percentileLims = [1 99];  % remove and interpolate tracking outside this percentile range
s.plotPredictors = true;

% initializations
% if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end % reassign settings passed in varargin
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'frameTimeStamps', 'frameTimeStampsWisk', 'wheelPositions', 'wheelTimes', ...
    'bodyAngles', 'whiskerAngle', 'pixelsPerM', 'wheelCenter', 'wheelRadius', ...
    'lickTimes');
locationsRun = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw.csv'));
locationsWisk = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw_wisk.csv'));
pawNames = {'paw1LH', 'paw2LF', 'paw3RF', 'paw4RH'};
dimensions = {'x', 'y', 'z'};
fps = 1/nanmedian(diff(frameTimeStamps));
scoreThresh = getScoreThresh(session, 'trackedFeaturesRaw_metadata.mat');  % scoreThresh depends on whether deeplabcut (old version) or deepposekit was used


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
pawXYZ = getPawXYZ(session);
for i = 1:length(pawNames)
    for j = 1:3
        pos = pawXYZ.(pawNames{i})(:,j);
        pos = interp1(frameTimeStamps, pos, t);
        vel = getVelocity(pos, s.velTime, fps);
        cntPredictors.([pawNames{i} '_' dimensions{j}]) = pos(:);
        cntPredictors.([pawNames{i} '_' dimensions{j} '_vel']) = vel(:);
    end
end

% body angle
cntPredictors.bodyAngle = interp1(frameTimeStamps, bodyAngles, t)';

% whisker angle
cntPredictors.whiskerAngle = interp1(frameTimeStamps, whiskerAngle, t)';

% butt height
[buttHeight, buttConfidence] = deal(locationsRun.tailBase_top_1, locationsRun.tailBase_top_2);
buttHeight(buttConfidence<scoreThresh) = nan;
buttHeight = fillmissing(buttHeight, 'pchip');
buttHeight = ((wheelCenter(2)-wheelRadius) - buttHeight) / pixelsPerM; % flip z and set s.t. top of wheel is zero
buttHeight = interp1(frameTimeStamps, buttHeight, t);
cntPredictors.buttHeight = buttHeight(:);

% jaw (first pc projection)
[jawX, jawZ, jawConfidence] = deal(locationsWisk.jaw, locationsWisk.jaw_1, locationsWisk.jaw_2);
jaw = pcProject([jawX, jawZ], jawConfidence);
jaw = (jaw-nanmean(jaw)) / pixelsPerM;
jaw = interp1(frameTimeStampsWisk, jaw, t);
cntPredictors.jaw = jaw(:);

% ear (first pc projection)
[earX, earZ, earConfidence] = deal(locationsRun.ear, locationsRun.ear_1, locationsRun.ear_2);
ear = pcProject([earX, earZ], earConfidence);
ear = (ear-nanmean(ear)) / pixelsPerM;
ear = interp1(frameTimeStamps, ear, t);
cntPredictors.ear = ear(:);

% nose (first pc projection)
[noseX, noseZ, noseConfidence] = deal(locationsWisk.nose, locationsWisk.nose_1, locationsWisk.nose_2);
nose = pcProject([noseX, noseZ], noseConfidence);
nose = (nose-nanmean(nose)) / pixelsPerM;
nose = interp1(frameTimeStampsWisk, nose, t);
cntPredictors.nose = nose(:);

% satiation
lickBins = histcounts(lickTimes, [t (t(end)+s.dt)]);
satiation = cumsum(lickBins);
satiation = satiation / max(satiation);
cntPredictors.satiation = satiation(:);



% -------
% LOGICAL
% -------






% make nice little plot
if s.plotPredictors
    dz = 8;  % (standard deviation) vertical separation
    xLims = [992 1042];
    
    figure('color', 'white', 'position', [125.00 93.00 1666.00 886.00]);
    names = cntPredictors.Properties.VariableNames(2:end);
    plot(t, zscore(table2array(cntPredictors(:,2:end)))+ [1:width(cntPredictors)-1]*dz);
    set(gca, 'xlim', xLims, 'visible', 'off')
    for i = 1:length(names)
        text(xLims(1), i*dz, names{i}, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'Interpreter', 'none')
    end
end



function projection = pcProject(sig, confidence)
    % project signal onto first principal component
    % signal is (time X dimension) matrix
    
    sig(confidence<scoreThresh,:) = nan;
    pcs = pca(sig);
    projection = sig * pcs(:,1);  % project onto first pc
    lims = prctile(projection, s.percentileLims);
    projection(projection<lims(1) | projection>lims(2)) = nan;
    projection = fillmissing(projection, 'pchip');
end

end




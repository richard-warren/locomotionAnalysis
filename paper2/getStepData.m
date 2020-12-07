function stepData = getStepData(session, varargin)
% for each paw, gets inds of each step and metadata associated with each
% step, e.g. wheel vel, body angle, isStepOverObs, paw height, etc... //
% this table will subsequently be used to compute step PSTHs that are
% binned by various step characteristics


% settings
s.tlims = [.1 .5];      % (s) exclude steps shorter than tlims(1) or longer than tlims(2)
s.velTime = .01;        % (s) compute wheel velocity over this time window
s.outputFileName = '';  % if provided saves stepData to disk


% load session data
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'frameTimeStamps', 'wheelPositions', 'wheelTimes', 'wheelCenter', 'wheelRadius', ...
    'pixelsPerM', 'bodyAngles', 'obsPixPositions', 'obsOnTimes', 'obsOffTimes');
wheelFs = 1/median(diff(wheelTimes));
fps = 1/median(diff(frameTimeStamps));  % video frames per second
vel = getVelocity(wheelPositions, s.velTime, wheelFs);

% get stance bins
pawNames = {'paw1LH', 'paw2LF', 'paw3RF', 'paw4RH'};
[pawXYZ, pawXYZ_pixels] = getPawXYZ(session);
xzLocations = nan(length(frameTimeStamps), 2, 4);  % (time X xz X paw)
for i = 1:length(pawNames); xzLocations(:,:,i) = pawXYZ_pixels.(pawNames{i})(:,[1,3]); end  % format required by getStanceBins()
stanceBins = getStanceBins(frameTimeStamps, xzLocations, wheelPositions, wheelTimes, wheelCenter, wheelRadius, fps, pixelsPerM);

% compute 'unheadfixed' x coords
wheelPosInterp = interp1(wheelTimes, wheelPositions, frameTimeStamps);  % wheel positions at frame times
xunheadfixed = nan(length(frameTimeStamps), 4);
for i = 1:4; xunheadfixed(:,i) = pawXYZ{:,i+1}(:,1); end
xunheadfixed = xunheadfixed + wheelPosInterp;


% initialize tables (one per paw)
stepData = cell(1,4);
for i = 1:4
    % get step inds
    temp = find(diff(stanceBins(:,i))==1);  % inds right before a stances start
    startInds = temp(1:end-1)+1;            % inds for stance starts
    endInds = temp(2:end);                  % ind before start of next stance
    
    % remove outlier steps
    dts = (endInds-startInds) / fps;
    bins = dts>=s.tlims(1) & dts<=s.tlims(2);
    
    % initialize table
    inds = [startInds(bins), endInds(bins)];
    times = frameTimeStamps(inds);
    n = sum(bins);
    
    vars =  {'times', 'inds', 'vel',    'bodyAngle', 'isStepOver', 'length', 'height'};
    inits = {times,   inds,   nan(n,1), nan(n,1),    false(n,1),   nan(n,1), nan(n,1)};
    stepData{i} = table(inits{:}, 'VariableNames', vars);
    
end

% compute step metadata
for i = 1:4
    
    times = stepData{i}.times;
    inds = stepData{i}.inds;
    
    % velocity
    stepData{i}.vel = interp1(wheelTimes, vel, times(:,1));
    
    % body angle
    stepData{i}.bodyAngle = interp1(frameTimeStamps, bodyAngles, times(:,1));
    
    % length
    stepData{i}.length = diff(xunheadfixed(inds), 1, 2);
    
    % height
    zvals = pawXYZ{:,i+1}(:,3);
    for j = 1:size(times,1)
        stepData{i}.height(j) = max(zvals(inds(j,1):inds(j,2)));
    end
    
    % isStepOver
    obsEndPos = interp1(frameTimeStamps, obsPixPositions, times(:,2));             % x positions of obs (in pixels) at the end of every step
    pawEndPos = interp1(frameTimeStamps, pawXYZ_pixels{:, i+1}(:,1), times(:,2));  % x positions of paw (in pixels) at the end of every step
    for j = 1:length(obsOnTimes)
        stepOverInd = find(times(:,2) > obsOnTimes(j) & ...
                           times(:,2) < (obsOffTimes(j)+1) & ...  % add one second buffer
                           pawEndPos  > obsEndPos, 1, 'first');
        stepData{i}.isStepOver(stepOverInd) = true;
    end
end

% save
if ~isempty(s.outputFileName)
    settings = s;
    save(s.outputFileName, 'stepData', 'settings');
end






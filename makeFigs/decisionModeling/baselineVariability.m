%% load experiment data
fprintf('loading... '); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'baseline_data.mat'), 'data'); disp('baseline data loaded!')

% initializations
mice = {data.data.mouse};
global_config;


%% ----------
% PLOT THINGS
%  ----------

%% distance and time to contact

% settings
trialSmps = 100;

% initializations
[distances, times] = deal(cell(1,length(mice)));

for i = 1:length(mice)
    
    fprintf('%s: collecting data...\n', mice{i})
    sessions = {data.data(i).sessions.session};  % sessions for mouse
    
    [distances{i}, times{i}] = deal([]);
    for j = 1:length(sessions)
        
        % load session data
        load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{j}, 'kinData.mat'), 'kinData')
        load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{j}, 'runAnalyzed.mat'), 'frameTimeStamps', 'isLightOn');
        secondsPerFrame = nanmedian(diff(frameTimeStamps)); % seconds per frame
        
        % set conditionals here
        modPawOnlySwing = false(size(kinData)); % only trials where mod paw is in swing and non mod paw is in stance
        modPawOnlySwing([kinData.isTrialAnalyzed]) = ([kinData.isLeftSwingAtContact] + [kinData.isRightSwingAtContact]) == 1;
        bins = [kinData.isTrialAnalyzed] & modPawOnlySwing & ~isLightOn(:)';
        
        for k = find(bins)
            
            % get distance of leading paw at contact
            contactInd = find(frameTimeStamps(kinData(k).trialInds) >= kinData(k).wiskContactTimes, 1, 'first'); % ind in trial at which contact occurs
            distances{i}(end+1) = abs(max([kinData(k).locations(contactInd,1,:)]))*1000;

            trialX = max(kinData(k).locations(contactInd-trialSmps+1:contactInd,1,:), [], 3);
            linFit = polyfit(trialX', 1:trialSmps, 1);
            predictedAtObsInd = polyval(linFit, 0);
            times{i}(end+1) = abs((predictedAtObsInd-trialSmps) * secondsPerFrame)*1000; % frame until contact * (seconds/frame)
            
        end
        clear kinData frameTimeStamps isLightOn
    end
end
disp('all done!')

%% PLOT DISTANCE AND TIME TO CONTACT

% settings
xLims = [15 50];
yLims = [0 150];
scatSize = 4;
scatAlpha = .5;
mouseColors = true;
scatPlotSize = .7;
border = .15;


% initializations
d = abs(cat(2, distances{:}));
t = abs(cat(2, times{:}));
mouseIds = repelem(1:length(mice), cellfun(@length, distances));
medDistances = cellfun(@nanmedian, distances);
medTimes = cellfun(@nanmedian, times);

% plot that shit
figure('Color', 'white', 'Position', [2000 400 450 350], 'MenuBar', 'none');
scatterHistoRick(d,t, ...
    {'groupId', mouseIds, 'colors', 'jet', ...
    'xlabel', 'distance to contact (mm)', 'ylabel', 'time to contact (ms)', ...
    'xLims', xLims, 'yLims', yLims, 'showCrossHairs', true, 'scatSize', scatSize, 'scatAlpha', scatAlpha});


% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'baseline_distanceTimeToContact');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');

%% SCATTER VEL AND BODY ANGLE AT CONTACT

flat = flattenData(data, {'mouse', 'session', 'trial', 'velAtWiskContact', 'angleAtWiskContact', 'isLightOn', 'modPawOnlySwing'});
flat = flat(~[flat.isLightOn]);
flat = flat([flat.modPawOnlySwing]);
[~,~,mouseIds] = unique({flat.mouse});

% settings
xLims = [0 1];
yLims = [-30 30];

% initializations
figure('Color', 'white', 'Position', [2000 400 450 350], 'MenuBar', 'none');
scatterHistoRick([flat.velAtWiskContact], [flat.angleAtWiskContact], ...
    {'groupId', mouseIds, 'colors', 'jet', ...
    'xlabel', 'velocity (m/s)', 'ylabel', ['body angle (' char(176) ')'], ...
    'xLims', xLims, 'yLims', yLims, 'showCrossHairs', true, 'scatSize', scatSize, 'scatAlpha', scatAlpha});

% save
file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'matlabFigs', ...
        'baseline_phaseVelVariability');
fprintf('writing %s to disk...\n', file)
saveas(gcf, file, 'svg');












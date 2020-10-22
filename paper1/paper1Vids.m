%% inits

paper1_config;
movieDir = fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'movies');


%% make wide view setup example view (shows obs coming into view from far off on the right side of the screen)
% (movie 1)

% settings
session = '180913_003';
trials = [7 10 30];  % [7 10 30 50]
speed = .25;

makeSetupExampleVid(fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'movies', 'setupExample.mp4'), ...
    session, trials, 'includeWiskCam', true, 'obsColor', obsColor, 'obsOffsetX', 1, 'obsOffsetZ', 6, ...
    'playBackSpeed', speed, 'text', sprintf('%i%% speed', speed*100), ...
    'textArgs', {'FontSize', 20, 'Font', 'Arial'});

%% tracking example vid
% (movie 1)

% settings
session = '180703_002';
trials = [74 83 96];
speed = .1;

makeVid(fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'movies', 'trackingExample'), ...
    session, 'trials', trials, 'playBackSpeed', speed, 'obsPosRange', [-.15 .15], ...
    'showTracking', true, 'text', sprintf('%i%% speed', speed*100), 'colors', trackingColors, ...
    'textArgs', {'FontSize', 12, 'Font', 'Arial'}, 'insertTrialInfo', false)


% % uncomment to make options
% numTrials = 10;
% sessions = getAllExperimentSessions('experiments', 'baselineNotes'); sessions = sessions.session;
% for i = 1:length(sessions)
%     try
%         % find fastest numTrials trials in session
%         load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'runAnalyzed.mat'), 'obsOnTimes', 'obsOffTimes')
%         [~, sortInds] = sort(obsOffTimes - obsOnTimes);
%         trials = sort(sortInds(1:numTrials)');
%         
%         makeVid(fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'movies', 'movie_options', 'trackingExample', ['trackingExample_' sessions{i}]), ...
%             sessions{i}, 'trials', trials, 'playBackSpeed', speed, 'obsPosRange', [-.1 .1], ...
%             'showTracking', true, 'text', sprintf('%i%% speed', speed*100), 'colors', trackingColors, ...
%             'textArgs', {'FontSize', 12, 'Font', 'Arial'}, 'insertTrialInfo', true)
%     catch
%         fprintf('PROBLEM WITH SESSION %s\n', sessions{i})
%     end
% end


%% unheadfixed tracking
% (movie 1)

% settings
session = '180703_002';
trials = [80 96];
speed = .15;
dropFrames = 2;

makeVidUnheadfixed(fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'movies', 'trackingUnheadfixed'), ...
    session, 'trials', trials, 'playBackSpeed', speed, 'showTracking', true, 'dropFrames', dropFrames, ...
    'colors', cat(1, stepColors, ones(2,3)*.5), 'text', sprintf('%i%% speed', speed*100), ...
    'textArgs', {'FontSize', 20, 'FontName', 'Arial'}, 'insertTrialInfo', false)


%% mice clear obstacles at high speeds
% (movie 2)

% settings
session = '180715_004';
trials = [8 25];
slowDownSpeed = .1;  % playback speed right after whisker contact
insertTrialInfo = false;

% make real-time vid
makeVid(fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'movies', 'highSpeedObstacleClearance'), ...
    session, 'trials', trials(1), 'playBackSpeed', 1.0, 'rewardTrials', true, 'rewardWindow', [5 -1], 'textArgs', {'FontSize', 12, 'Font', 'Arial'}, ...
    'visible', 'on', 'showTracking', true, 'text', '100% speed', 'dropFrames', 5, 'colors', trackingColors, 'insertTrialInfo', insertTrialInfo)

% make real-time vid that slows at whisker contact
makeVid(fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'movies', 'highSpeedObstacleClearance_slowDown'), ...
    session, 'trials', trials, 'playBackSpeed', 1.0, 'rewardTrials', true, 'rewardWindow', [5 -1], 'textArgs', {'FontSize', 12, 'Font', 'Arial'}, ...
    'visible', 'on', 'showTracking', true, 'text', sprintf('100%% speed, %i%% speed during obstacle clearance', slowDownSpeed*100), ...
    'dropFrames', 5, 'speedNearContact', slowDownSpeed, 'colors', trackingColors, 'insertTrialInfo', insertTrialInfo, 'contactWindow', [-.02 .2])

% makeVidUnheadfixed(fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'movies', 'highSpeedObstacleClearance_unheadfixed'), ...
%     session, 'trials', trials, 'playBackSpeed', 1.0, 'showTracking', false, 'dropFrames', 5, 'text', '100% speed', ...
%     'colors', trackingColors, 'insertTrialInfo', insertTrialInfo)


%% sensory dependence
% (movie 3)

% settings
speed = .15;
vidArgs = {'playBackSpeed', speed, 'showTracking', false, 'text', sprintf('%i%% speed', speed*100), ...
    'dropFrames', 2, 'insertTrialInfo', false, 'obsPosRange', [-.15 .1], 'xLims', [-.2 .1], 'textArgs', {'FontSize', 14, 'FontName', 'Arial'}};

noWiskSession = '180805_004';  % mtc5
wiskSession = '180803_004';    % mtc5
trials_WV = [37 52 100];
trials_W  = [144 39 45];
trials_V  = [196 32 64];
trials_neither  = [19 29 39];

% WV
makeVidUnheadfixed(fullfile(movieDir, 'sensoryDependence_WV'), wiskSession,   'trials', trials_WV, vidArgs{:});
makeVidUnheadfixed(fullfile(movieDir, 'sensoryDependence_W'),  wiskSession,   'trials', trials_W, vidArgs{:});
makeVidUnheadfixed(fullfile(movieDir, 'sensoryDependence_V'),  noWiskSession, 'trials', trials_V, vidArgs{:});
makeVidUnheadfixed(fullfile(movieDir, 'sensoryDependence_-'),  noWiskSession, 'trials', trials_neither, vidArgs{:});


% uncomment to make options
% randomSeed = 0;
% numTrials = 20;
% miceToExclude = {'sen1', 'cmu6'};
% 
% sessions = getAllExperimentSessions('experiments', 'sensoryDependenceNotes', 'extraColumns', {'whiskers'});
% sessions = sessions(~ismember(sessions.mouse, miceToExclude), :);
% mice = unique(sessions.mouse);
% 
% 
% for i = 1:length(mice)
%     
%     % select sessions and load data
%     wiskSessions   = sessions(strcmp(sessions.mouse, mice{i}) & strcmp(sessions.whiskers, 'full'), :).session;
%     noWiskSessions = sessions(strcmp(sessions.mouse, mice{i}) & strcmp(sessions.whiskers, 'none'), :).session;
%     wiskSession = wiskSessions{1};
%     noWiskSession = noWiskSessions{1};
%     wiskData   = load(fullfile(getenv('OBSDATADIR'), 'sessions', wiskSession, 'runAnalyzed.mat'), 'isLightOn', 'isWheelBreak', 'obsHeights');
%     noWiskData = load(fullfile(getenv('OBSDATADIR'), 'sessions', noWiskSession, 'runAnalyzed.mat'), 'isLightOn', 'isWheelBreak', 'obsHeights');
%     
%     % WV
%     if ~exist(fullfile(getenv('OBSDATADIR'), 'sessions', wiskSession, 'run.mp4'), 'file'); concatTopBotVids(wiskSession); end  % temp
%     trials = find(wiskData.isLightOn & ~wiskData.isWheelBreak & ~isnan(wiskData.obsHeights)');
%     rng(randomSeed); trials = sort(datasample(trials, numTrials))';
%     makeVidUnheadfixed(fullfile(movieDir, 'movie_options', 'sensoryDependence', [mice{i} '_WV']), wiskSession, 'trials', trials, vidArgs{:}, 'insertTrialInfo', true);
%     
%     % W
%     trials = find(~wiskData.isLightOn & ~wiskData.isWheelBreak & ~isnan(wiskData.obsHeights)');
%     rng(randomSeed); trials = sort(datasample(trials, numTrials))';
%     makeVidUnheadfixed(fullfile(movieDir, 'movie_options', 'sensoryDependence', [mice{i} '_W']), wiskSession, 'trials', trials, vidArgs{:}, 'insertTrialInfo', true);
%     
%     % V
%     if ~exist(fullfile(getenv('OBSDATADIR'), 'sessions', noWiskSession, 'run.mp4'), 'file'); concatTopBotVids(noWiskSession); end  % temp
%     trials = find(noWiskData.isLightOn & ~noWiskData.isWheelBreak & ~isnan(noWiskData.obsHeights)');
%     rng(randomSeed); trials = sort(datasample(trials, numTrials))';
%     makeVidUnheadfixed(fullfile(movieDir, 'movie_options', 'sensoryDependence', [mice{i} '_V']), noWiskSession, 'trials', trials, vidArgs{:}, 'insertTrialInfo', true);
%     
%     % -
%     trials = find(~noWiskData.isLightOn & ~noWiskData.isWheelBreak & ~isnan(noWiskData.obsHeights)');
%     rng(randomSeed); trials = sort(datasample(trials, numTrials))';
%     makeVidUnheadfixed(fullfile(movieDir, 'movie_options', 'sensoryDependence', [mice{i} '_-']), noWiskSession, 'trials', trials, vidArgs{:}, 'insertTrialInfo', true);
% end




%% decision kinematics
% (movie 4)

% settings
session = '181216_003';  % 181216_003, 180730_001
speed = .15;

makeDecisionVid(fullfile(movieDir, 'decision'), session, ...
    'speed', speed, 'text', sprintf('%i%% speed', speed*100), 'dropFrames', 1, 'xLims', [-.18 .1], ...
    'textArgs', {'FontSize', 14, 'FontName', 'Arial'}, 'colors', [decisionColors; .8 .8 .8]);


% uncomment to make options
% sessions = getAllExperimentSessions('experiments', 'baselineNotes'); sessions = sessions.session;
% for i = 1:length(sessions)
%     try
%         makeDecisionVid(fullfile(movieDir, 'movie_options', 'decision', ['decision_' sessions{i}]), ...
%             sessions{i}, 'speed', speed, 'text', sprintf('%i%% speed', speed*100), 'dropFrames', 1, ...
%             'xLims', [-.2 .1], 'textArgs', {'FontSize', 14, 'FontName', 'Arial'}, 'colors', decisionColors);
%     catch
%         fprintf('PROBLEM WITH SESSION %s\n', sessions{i})
%     end
% end


%% barrel cortex lesions
% (movie 5)

% settings
speed = .15;
vidArgs = {'playBackSpeed', speed, 'showTracking', false, 'text', sprintf('%i%% speed', speed*100), ...
    'dropFrames', 2, 'insertTrialInfo', false, 'obsPosRange', [-.15 .1], 'xLims', [-.15 .1], 'textArgs', {'FontSize', 14, 'FontName', 'Arial'}};
preSession  = '191106_001';  % sen13
postSession = '191107_001';  % sen13
preTrials  = [40 157 234];
postTrials = [44 86 128];

makeVidUnheadfixed(fullfile(movieDir, 'senLesion_pre'),  preSession,  'trials', preTrials,  vidArgs{:});
makeVidUnheadfixed(fullfile(movieDir, 'senLesion_post'), postSession, 'trials', postTrials, vidArgs{:});


% uncomment to make options
% sessions = getAllExperimentSessions('experiments', 'senLesionNotes', 'extraColumns', {'condition'});
% mice = unique(sessions.mouse);
% numTrials = 10;
% for i = 1:length(mice)
%     preInd  = find(strcmp(sessions.mouse, mice{i}) & ismember(sessions.condition, {'pre', 'postIpsi'}), 1, 'last');        % sessions ind for final 'pre' session
%     postInd = find(strcmp(sessions.mouse, mice{i}) & ismember(sessions.condition, {'postContra', 'postBi'}), 1, 'first');  % sessions ind for first 'post' session
%     if (postInd-preInd) ~= 1; disp('WARNING! Something fishy about session determination!'); end
% 
%     makeVidUnheadfixed(fullfile(movieDir, 'movie_options', 'barrel_cortex', [mice{i} '_pre']),  sessions.session{preInd},  'numTrials', numTrials, vidArgs{:}, 'insertTrialInfo', true);
%     makeVidUnheadfixed(fullfile(movieDir, 'movie_options', 'barrel_cortex', [mice{i} '_post']), sessions.session{postInd}, 'numTrials', numTrials, vidArgs{:}, 'insertTrialInfo', true);
% end

%% motor cortex lesions
% (movie 6)

% settings
speed = .15;
vidArgs = {'playBackSpeed', speed, 'showTracking', false, 'text', sprintf('%i%% speed', speed*100), ...
    'dropFrames', 2, 'insertTrialInfo', false, 'obsPosRange', [-.15 .1], 'xLims', [-.15 .1], 'textArgs', {'FontSize', 14, 'FontName', 'Arial'}};
preSession  = '180827_001';  % mtc1
postSession = '180828_001';  % mtc1
preTrials  = [132 155 165];
postTrials = [91 118 207];

makeVidUnheadfixed(fullfile(movieDir, 'mtc_pre'),  preSession,  'trials', preTrials,  vidArgs{:});
makeVidUnheadfixed(fullfile(movieDir, 'mtc_post'), postSession, 'trials', postTrials, vidArgs{:});


% uncomment to make options
% sessions = getAllExperimentSessions('experiments', 'mtcLesionNotes', 'extraColumns', {'condition', 'brainRegion'});
% mice = unique(sessions.mouse);
% numTrials = 10;
% for i = 1:length(mice)
%     preInd  = find(strcmp(sessions.mouse, mice{i}) & strcmp(sessions.condition, 'pre'), 1, 'last');        % sessions ind for final 'pre' session
%     postInd = find(strcmp(sessions.mouse, mice{i}) & strcmp(sessions.condition, 'post'), 1, 'first');  % sessions ind for first 'post' session
%     if (postInd-preInd) ~= 1; disp('WARNING! Something fishy about session determination!'); end
%     makeVidUnheadfixed(fullfile(movieDir, 'movie_options', 'motor_cortex', 'lesions', [mice{i} '_' sessions.session{preInd} '_pre']),   sessions.session{preInd},  'numTrials', numTrials, vidArgs{:}, 'insertTrialInfo', true);
%     makeVidUnheadfixed(fullfile(movieDir, 'movie_options', 'motor_cortex', 'lesions', [mice{i} '_' sessions.session{postInd} '_post']), sessions.session{postInd}, 'numTrials', numTrials, vidArgs{:}, 'insertTrialInfo', true);
% end

%% motor cortex muscimol
% (movie 6)

% settings
speed = .15;
vidArgs = {'playBackSpeed', speed, 'showTracking', false, 'text', sprintf('%i%% speed', speed*100), ...
    'dropFrames', 2, 'insertTrialInfo', false, 'obsPosRange', [-.15 .1], 'xLims', [-.15 .1], 'textArgs', {'FontSize', 14, 'FontName', 'Arial'}};
salSession = '180716_003';  % mtc1
musSession = '180717_005';  % mtc1
salTrials = [49 71 84];
musTrials = [7 18 84];

makeVidUnheadfixed(fullfile(movieDir, 'mtc_saline'),   salSession, 'trials', salTrials, vidArgs{:});
makeVidUnheadfixed(fullfile(movieDir, 'mtc_muscimol'), musSession, 'trials', musTrials, vidArgs{:});



% uncomment to make options
% sessions = getAllExperimentSessions('experiments', 'muscimolNotes', 'extraColumns', {'condition', 'brainRegion'});
% mice = unique(sessions.mouse);
% numTrials = 10;
% for i = 1:length(mice)
%     salInd  = find(strcmp(sessions.mouse, mice{i}) & strcmp(sessions.condition, 'saline'), 1, 'first');     % sessions ind for first 'saline' session
%     musInd = find(strcmp(sessions.mouse, mice{i}) & strcmp(sessions.condition, 'muscimol'), 1, 'first');    % sessions ind for first 'muscimol' session
%     makeVidUnheadfixed(fullfile(movieDir, 'movie_options', 'motor_cortex', 'muscimol', [mice{i} '_' sessions.session{salInd} '_saline']),   sessions.session{salInd}, 'numTrials', numTrials, vidArgs{:}, 'insertTrialInfo', true);
%     makeVidUnheadfixed(fullfile(movieDir, 'movie_options', 'motor_cortex', 'muscimol', [mice{i} '_' sessions.session{musInd} '_muscimol']), sessions.session{musInd}, 'numTrials', numTrials, vidArgs{:}, 'insertTrialInfo', true);
% end






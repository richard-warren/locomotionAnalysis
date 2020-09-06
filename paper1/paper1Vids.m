%% inits

paper1_config;


%% make wide view setup example view (shows obs coming into view from far off on the right side of the screen)
% (movie 1)

% settings
session = '180913_003';
trials = [7 10 30 50];
speed = .15;

makeSetupExampleVid(fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'movies', 'setupExample.mp4'), ...
    session, trials, 'includeWiskCam', true, 'obsColor', obsColor, 'obsOffsetX', -20, 'obsOffsetZ', 6, ...
    'playBackSpeed', speed, 'text', sprintf('%.2fx', speed));

%% tracking example vid
% (movie 2)

% settings
session = '180715_004';
trials = 20:30;
speed = .1;

makeVid(fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'movies', 'trackingExample'), ...
    session, 'trials', trials, 'playBackSpeed', speed, 'obsPosRange', [-.1 .1], ...
    'showTracking', true, 'text', sprintf('%.2fx', speed))

%% mice clear obstacles at high speeds
% (movie 3)

% settings
session = '180715_004';
numTrials = 5;
slowDownSpeed = .1;  % playback speed right after whisker contact

% find fastest numTrials trials in session
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), 'rewardTimes')
dts = diff(rewardTimes);
[~, sortInds] = sort(dts);
trials = sort(sortInds(1:numTrials)'+1);

% make real-time vid
makeVid(fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'movies', 'highSpeedObstacleClearance'), ...
    session, 'trials', trials, 'playBackSpeed', 1.0, 'rewardTrials', true, 'rewardWindow', [5 -1], ...
    'visible', false, 'showTracking', true, 'text', sprintf('%.2fx', 1.0), 'dropFrames', 5)

% make real-time vid that slows at whisker contact
makeVid(fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'movies', 'highSpeedObstacleClearance_slowDown'), ...
    session, 'trials', trials, 'playBackSpeed', 1.0, 'rewardTrials', true, 'rewardWindow', [5 -1], ...
    'visible', false, 'showTracking', true, 'text', sprintf('%.2fx, %.2fx at whisker contact', 1.0, slowDownSpeed), ...
    'dropFrames', 5, 'speedNearContact', slowDownSpeed)

makeVidUnheadfixed(fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'movies', 'highSpeedObstacleClearance_unheadfixed'), ...
    session, 'trials', trials, 'playBackSpeed', 1.0, 'showTracking', false, 'dropFrames', 5, 'text', sprintf('%.2fx', 1.0))

%% unheadfixed tracking
% (movie 4)

% settings
session = '180715_004';
trials = [7 10 30 50];
speed = .15;
dropFrames = 2;

makeVidUnheadfixed(fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'movies', 'trackingUnheadfixed'), ...
    session, 'trials', trials, 'playBackSpeed', speed, 'showTracking', true, 'dropFrames', dropFrames, ...
    'colors', cat(1, stepColors, ones(2,3)*.5), 'text', sprintf('%.2fx', speed))


%% decision kinematics
% (movie 6)

% settings
session = '180630_000';
speed = .15;

% makeDecisionVid(fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'movies', 'decision'), ...
%     session, 'speed', speed, 'text', sprintf('%.2fx', speed), 'dropFrames', 1, 'xLims', [-.25 .1]);


% uncomment to make options
sessions = getAllExperimentSessions('experiments', 'baselineNotes'); sessions = sessions.session;
for i = 1:length(sessions)
    if ~exist(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'run.mp4'), 'file'); concatTopBotVids(sessions{i}); end
    try
        makeDecisionVid(fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'movies', 'movie_options', ['decision_' sessions{i}]), ...
            sessions{i}, 'speed', speed, 'text', sprintf('%.2fx', speed), 'dropFrames', 1, 'xLims', [-.25 .1]);
    catch
        fprintf('PROBLEM WITH SESSION %s\n', sessions{i})
    end
end












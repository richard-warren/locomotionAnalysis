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

%% tracking example vid (movie 2)

% settings
session = '180715_004';
trials = 20:30;
speed = .1;

makeVid(fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'movies', 'trackingExample'), ...
    session, 'trials', trials, 'playBackSpeed', speed, 'obsPosRange', [-.1 .1], ...
    'showTracking', true, 'text', sprintf('%.2fx', speed))

%% mice clear obstacles at high speeds (movie 3)

% settings
session = '180715_004';
trials = [];
speed = 1.0;

makeVid(fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'movies', 'highSpeedObstacleClearance'), ...
    session, 'trials', trials, 'playBackSpeed', speed, 'obsPosRange', [-.1 .1], ...
    'showTracking', true, 'text', sprintf('%.2fx', speed))






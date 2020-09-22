%% inits

stepSaturation = .9;
stepColors = hsv2rgb([.65 stepSaturation 1;
                      .55 stepSaturation 1;
                      .02 stepSaturation 1;
                      .12 stepSaturation 1]);  % LF, TF, LH, TH
trackingColors = cat(1, stepColors, ones(2,3)*.5);  % colors for tracking (first four rows are paw, final two are tail)


%% tracking examples, both fast and slowed down
% (movie 1)

% settings
close all
session = '191113_000';
% trials = [74 83 96];
trialNum = 5;
speed = .1;
obsLims = [-1.1 -.06];
contrast = [.1 .75];

% pick fastest trialNum trials
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), 'obsOnTimes', 'obsOffTimes')
dts = obsOffTimes - obsOnTimes;
[~, trials] = sort(dts);
trials = sort(trials(1:trialNum))';

for includeWiskCam = [false true]
    if includeWiskCam; suffix=''; else; suffix='_noWiskCam'; end

    makeVid(fullfile(getenv('OBSDATADIR'), 'papers', 'wheel_paper', ['tracking_fast' suffix]), ...
        session, 'trials', trials, 'playBackSpeed', 1, 'obsPosRange', obsLims, ...
        'showTracking', true, 'text', '100% speed', 'colors', trackingColors, 'contrastLims', contrast, ...
        'textArgs', {'FontSize', 12, 'Font', 'Arial'}, 'insertTrialInfo', false, 'dropFrames', 5, 'includeWiskCam', includeWiskCam)

    makeVid(fullfile(getenv('OBSDATADIR'), 'papers', 'wheel_paper', ['tracking_slow' suffix]), ...
        session, 'trials', trials(1), 'playBackSpeed', speed, 'obsPosRange', obsLims, ...
        'showTracking', true, 'text', sprintf('%i%% speed', speed*100), 'colors', trackingColors, 'contrastLims', contrast, ...
        'textArgs', {'FontSize', 12, 'Font', 'Arial'}, 'insertTrialInfo', false, 'includeWiskCam', includeWiskCam)
end


%% unheadfixed tracking
% (movie 2)

% settings
session = '180703_002';
trials = [80 96];
speed = .15;
dropFrames = 2;
xLims = [-.6 -.1];

for includeWiskCam = [false true]
    if includeWiskCam; suffix=''; else; suffix='_noWiskCam'; end
    
    makeVidUnheadfixed(fullfile(getenv('OBSDATADIR'), 'papers', 'wheel_paper', ['unheadfixed' suffix]), ...
        session, 'trials', trials, 'playBackSpeed', speed, 'showTracking', true, 'dropFrames', dropFrames, ...
        'colors', cat(1, stepColors, ones(2,3)*.5), 'text', sprintf('%i%% speed', speed*100), ...
        'textArgs', {'FontSize', 20, 'FontName', 'Arial'}, 'insertTrialInfo', false, 'xLims', xLims, 'includeWiskCam', includeWiskCam)
end


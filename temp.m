session = '190205_001';

load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'run.mat'), 'whEncodA')
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), 'wheelPositions', 'wheelTimes')

wheelCircumference = 0.5984734005; % meters
encoderTicks = 1440; % how many encoder ticks per rotation

times = whEncodA.times;
diffs = diff(times);
vel = repmat((wheelCircumference/encoderTicks), 1, length(diffs)) ./ diffs';
velMl = getVelocity(wheelPositions, .02, 1000);


fprintf('mean: %.4f\nmatlab mean: %.4f\n', mean(vel), nanmean(abs(velMl)));

close all; figure; plot(times(1:end-1), vel); hold on; plot(wheelTimes, abs(velMl))
figure; histogram(vel); hold on; histogram(velMl);
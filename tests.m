%% test getFrameTimes with fake data

while true

    % !!! breaks when misses frame right before reward...
    
    % settings
    frames = 250000;
    missedFrames = 100;
    rewards = 50;
    beginningDrop = 250;  % how many frames to drop at the beginning

    % simulate data
    ttlTimes = 100:(1/250):(frames/250+100);
    rewardInds = linspace(1, length(ttlTimes), rewards+2);
    rewardInds = round(rewardInds(2:end-1));
    deltas = zeros(size(ttlTimes));
    deltas(rewardInds) = .5;
    ttlTimes = ttlTimes + cumsum(deltas);

    frameTimesRaw = ttlTimes + 1000;  % simulate clock offsets
    frameCounts = 200:(200+length(frameTimesRaw));

    % drop frames
    bins = true(size(frameTimesRaw));
    missedInds = sort(randsample(length(frameTimesRaw)-beginningDrop+1, missedFrames)+beginningDrop+1);
    bins(missedInds) = false;
    frameTimesRaw = frameTimesRaw(bins);
    frameCounts = frameCounts(bins);

    % get ground truth frame times
    frameTimesTrue = ttlTimes';
    frameTimesTrue = frameTimesTrue(bins);  % remove missing frames

    % remove first second of frames
    frameTimesRaw = frameTimesRaw(beginningDrop+1:end);
    frameCounts = frameCounts(beginningDrop+1:end);
    frameTimesTrue = frameTimesTrue(beginningDrop+1:end);

%     frameTimesTrue(1:find(diff(frameTimesTrue)>.1,1,'first')) = nan; % replacing first trial with nans

    % run algorithm
%     frameTimes = getFrameTimes2(ttlTimes', frameTimesRaw', frameCounts', 'test');
%     frameTimes = getFrameTimes4(ttlTimes', frameTimesRaw', frameCounts', 'test');
    frameTimes = getFrameTimes(ttlTimes', frameTimesRaw', frameCounts', 'test');

    % show where they differ
    close all; figure('position', [2402.00 276.00 560.00 420.00]);
    plot(frameTimesTrue, 'lineWidth', 2); hold on; plot(frameTimes)
    frameTimesTrue(isnan(frameTimesTrue)) = 0; frameTimes(isnan(frameTimes)) = 0;  % in matlab nan~=nan, so turn nans to zeros here for sake of comparison
    differInds = find(frameTimesTrue~=frameTimes);
    scatter(differInds, frameTimesTrue(differInds))
    fprintf('differing frames: %i\n', length(differInds))
    
    if ~isequal(frameTimes, frameTimesTrue); break; end
end


%% recompute time stamps

sessions = getAllExperimentSessions;

% !!! pick up with fixes for session 180807_000 (121)
for i = 84:length(sessions)
    
    fprintf('\n\n________________session #%i (%s)________________\n', i, experiments{i})
    load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'runAnalyzed.mat'), ...
        'frameTimeStamps', 'frameTimeStampsWisk')
    frameTimeStamps_old = frameTimeStamps;
    frameTimeStampsWisk_old = frameTimeStampsWisk;
    analyzeSession(sessions{i}, {'frameTimeStamps', 'frameTimeStampsWisk'})
    load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'runAnalyzed.mat'), ...
        'frameTimeStamps', 'frameTimeStampsWisk')
    
    
    differInds = find(frameTimeStamps ~= frameTimeStamps_old & ~isnan(frameTimeStamps_old));
    differIndsWisk = find(frameTimeStampsWisk ~= frameTimeStampsWisk_old & ~isnan(frameTimeStampsWisk_old));
    
    if ~isempty(differInds) || ~isempty(differIndsWisk)
        figure('Position', [1950.00 28.00 577.00 945.00], 'color', 'white', 'name', sessions{i});
    
        subplot(2,1,1); hold on
        plot(frameTimeStamps_old, 'linewidth', 2)
        plot(frameTimeStamps);

        scatter(differInds , frameTimeStamps(differInds))

        subplot(2,1,2); hold on
        plot(frameTimeStampsWisk_old, 'linewidth', 2)
        plot(frameTimeStampsWisk);
        scatter(differIndsWisk , frameTimeStampsWisk(differIndsWisk))
        pause(.1)
    end    
end


%%  test timeStampDecoderFLIR

% plots times decoded with this function to ensure they are all going in a straight line

sessions = getAllExperimentSessions;
times = cell(1, length(sessions));
close all; figure('Position', [2027.00 434.00 560.00 420.00]); hold on

for i = 1:length(sessions)
    disp(i/length(sessions))
    try
        camMetadata = dlmread(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'run.csv')); % columns: bonsai timestamps, point grey counter, point grey timestamps (uninterpretted)
        times{i} = timeStampDecoderFLIR(camMetadata(:,3));
        plot(times{i});
        pause(.01)
    end
end

%% test new rotary decoder

% plot old and new versions of decoded data

sessions = getAllExperimentSessions;
targetFs = 1000; % frequency that positional data will be resampled to
whEncoderSteps = 2880; % 720cpr * 4
wheelRad = 95.25; % mm
obEncoderSteps = 1000; % 250cpr * 4
obsRad = 96 / (2*pi); % radius of timing pulley driving belt of obstacles platform


for i = 1:length(sessions)
    try 
        load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'runAnalyzed.mat'), 'obsPositions', 'obsTimes'); % columns: bonsai timestamps, point grey counter, point grey timestamps (uninterpretted)
        load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'run.mat'), 'obEncodA', 'obEncodB');
        [posNew, timesNew] = rotaryDecoder(...
            obEncodA.times, obEncodA.level, obEncodB.times, obEncodB.level, obEncoderSteps, obsRad, targetFs, sessions{i});
        
        close all; figure('Position', [1921.00 1.00 1920.00 1003.00]); hold on
        plot(obsTimes, obsPositions, 'lineWidth', 2); hold on
        plot(timesNew, posNew, 'lineWidth', 1); hold on      
    end
end











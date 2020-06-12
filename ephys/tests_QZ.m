c = 1

fprintf('plotting cell %i... ', unit_ids(c))
figure('Name', sprintf('%s cell %i', session, unit_ids(c)), 'Visible', showFigures, ...
    'Color', 'white'); hold on
spkInds = allSpkInds{c};


% EXTRACT WAVEFORM ACROSS CHANNELS
spkIndsSub = spkInds(round(linspace(2, length(spkInds), spkNum))); % get subpopulation of spikes evenly spaced out
spkIndsSubAll = uint64(repmat(spkWindowInds,1,length(spkIndsSub)) + int64(repelem(spkIndsSub, length(spkWindowInds)))');

allWaveforms = getVoltage(data, 1:channelNum, spkIndsSubAll);
allWaveforms = reshape(allWaveforms, channelNum, length(spkWindowInds), []);
allWaveforms = permute(allWaveforms, [3 1 2]);
%     allWaveforms = allWaveforms - allWaveforms(:,:,1); % subtract beginning of trace from rest of trace (a hack of a high pass filter, lol)

% find best channel and get voltage for that channel
%     meanWaveform = squeeze(mean(allWaveforms,1));
%     [~, bestChannel] = max(peak2peak(meanWaveform,2));
channelData = getVoltage(data, bestChannels(c), 1:smps);
%     templates = readNPY(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, 'templates.npy'));
%     spike_templates = readNPY(fullfile(getenv('OBSDATADIR'), 'sessions', session, ephysFolder, 'spike_templates.npy'));



% PLOT SPIKE SHAPES ON PROBE
subplot(4,4,[1 5 9]); hold on
timeBins = discretize(double(spkIndsSub), timeBinNu409m);
colors = getColors(timeBinNum, cellColors(c,:));
sameShankInds = find(abs(xcoords - xcoords(bestChannels(c)))<50);

for j = sameShankInds
    for i = 1:timeBinNum
        firingRate = sum(timeBins==i) / (range(timeStamps)/timeBinNum);
        if firingRate>minFiringRate % don't plot average trace if rate of spikes in bin is too low, which happens when the unit is lost
            trace = squeeze(mean(allWaveforms(timeBins==i,j,:),1));
            plot(xcoords(j)*xSpacing + double(spkWindowInds), ...
                ycoords(j)*ySpacing + trace, ...
                'Color', colors(i,:), 'LineWidth', 2)
        end
    end
    if j==bestChannels(c); textColor='red'; else; textColor='black'; end
    text(double(xcoords(j)*xSpacing+spkWindowInds(1)), ...
        ycoords(j)*ySpacing, ...
        num2str(j), 'Color', textColor)
end
set(gca, 'visible', 'off')
    
%% inset figure

figure;
currentAxes(1) = subplot(4, 4, 1)
temp = [currentAxes(1).Position]
histogram(bodyAngles)

currentAxes(2) = axes('Position', [temp(1) + 0.115, temp(2) + 0.11, temp(3)*.2, temp(4)*.2 ])
box on
plot(linspace(1, 1000, length(bodyAngles)), zscore(bodyAngles));


%%

[bodyAnglesInterp, frameTimeStampsInterp, unitSpkRates, neuralTimesNew] = dataInterpolation(bodyAngles, frameTimeStamps, unit_spkRates, neuralTimes);


figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize'));

plotInd(1, 1) = 1;

plotInd(1, 2) = plotContinuousData(vel, unit_spkRates, 1, {'dataName', 'velocity (m/s)'});

plotInd(1, 3) = plotContinuousData(bodyAnglesInterp, unitSpkRates, plotInd(1, 2), {'dataName', 'bodyAngles'});


%%
figure;
[F,TF] = fillmissing(tailMid_bot(:, 2), 'spline');
x = 1:1000;
window = 1:1000;
plot(x, tailMid_bot(window, 2), 'b.-', x, F(window),'r.')
xlabel('frames');
ylabel('x position for tail mid point')
legend('Original Data','Filled Missing Data (spline)')


inds = intersect(find(~isnan(tailBase_bot(:, 1))), find(~isnan(tailBase_bot(:, 2))));






%% hilbert transform

% Generate sinusoid signal
dt = 0.001; % sampling interval [s]
t = dt:dt:1; % time axis [s]
f1 = 2.0; % freq of sinusoid [Hz]
phi0 = 0.0; % initial phase of sinusoid
d = sin(2.0*pi*t*f1 + phi0);

% Compute analytic signal
dA = hilbert(d); % The built -in " hilbert " function in MATLAB returns the analytic signal . We can plot the original signal , and the real and imaginary parts of the signal .

figure; clf() % clf = clear current figure window
subplot(4,1,1) % 4x1 plot , 1st plot
plot(d); ylabel('Data ');
subplot(4,1,2) % 4x1 plot , 2nd plot
plot(real(dA));
hold on;
plot(imag(dA), 'r');
hold off;
ylabel('Real ( blue ), Imag ( red)');
axis tight
% Note that the imaginary part of the analytic signal is the real part of the analytic signal shifted by 90 degrees (pi /2)
% Compute phase of analytic signal
phi = angle(dA);
% Compute amplitude envelope
amp = abs(dA);
% Plot results
subplot(4,1,3); plot(phi); ylabel('Angle '); axis tight
subplot(4,1,4); plot(amp); ylabel('Amplitude '); axis tight

%%
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'));
locationsTable = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw.csv')); % get raw tracking data
[locations, features] = fixTracking(locationsTable, frameTimeStamps, pixelsPerM);

pawsXLocations = locations(:, 1, 9:12);
pawNames = {'LH', 'LF', 'RF', 'RH'};


vel = getVelocity(wheelPositions, .02, 1000); % should have the same length as wheelTimes
wheelTimes = wheelTimes';

frameVel = nan(length(frameTimeStamps), 1);

temp = frameTimeStamps(:, 1);
[~,ind] = histc(temp,[-inf;(wheelTimes(1:end-1)+wheelTimes(2:end))/2;inf]);
frameTimeStamps(:, 2) = wheelTimes(ind(:, 1));
frameVel = vel(ind(:, 1))';

timeInds = find(frameVel>0.3);
selectedFrameTimeStamps = frameTimeStamps(timeInds);




figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on

paw{1} = pawsXLocations(timeInds, :, 2) - mean(pawsXLocations(timeInds, :, 2));
paw{2} = pawsXLocations(timeInds, :, 3) - mean(pawsXLocations(timeInds, :, 3));
% 
% paw{1} = medfilt1(paw{1}, 10);
% paw{2} = medfilt1(paw{2}, 10);

paw{1} = lowpass(paw{1}, 5, 250);
paw{2} = lowpass(paw{2}, 5, 250);

hPaw{1} = hilbert(squeeze(paw{1}));
hPaw{2} = hilbert(squeeze(paw{2}));

aPaw{1} = angle(hPaw{1}) + pi;
aPaw{2} = angle(hPaw{2}) + pi;


% aPaw{1} = medfilt1(aPaw{1}, 5);
% aPaw{2} = medfilt1(aPaw{2}, 5);

timePeriod = 240000:241000;

subplot(4, 1, 1);
x = aPaw{1}; y = aPaw{2}; 
plot(x(timePeriod), '-b', 'DisplayName', 'LF');hold on; axis tight;
plot(y(timePeriod), '-r', 'DisplayName', 'RF');
temp = aPaw{1} - aPaw{2}; 
phaseDiff = abs(mod(temp + pi, 2*pi) - pi);
plot(phaseDiff(timePeriod), '-k', 'DisplayName', 'Diff', 'LineWidth', 2);
legend({'LF', 'RF', 'Diff'})
xlabel('frameInds');
ylabel('phase (diff) after hilbert transformation');


subplot(4, 1, 2)
x = paw{1}; y = aPaw{1};
plot(x(timePeriod)/20 + 3, '-b'); axis tight; hold on
plot(y(timePeriod), '-r');
legend({'LF raw', 'LF hilbert transform'});
xlabel('frameInds');


subplot(4, 1, 3)
x = paw{2}; y = aPaw{2};
plot(x(timePeriod)/20 + 3, '-b'); axis tight; hold on
plot(y(timePeriod), '-r');
legend({'LF raw', 'LF hilbert transform'});
xlabel('frameInds');

subplot(4, 1, 4)
plot(zscore(phaseDiff), 'o');
xlabel('frameInds');
ylabel('zscores');



selectedInds = find(abs(zscore(phaseDiff)) > 4);
randomInds = datasample(selectedInds, 5, 'Replace', false);

for i = 1:length(randomInds)
    
    figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on

    
    timePeriod = randomInds(i) - 100 : randomInds(i) + 100;
    
    subplot(3, 1, 1);
    x = aPaw{1}; y = aPaw{2};
    plot(selectedFrameTimeStamps(timePeriod), x(timePeriod), '-b', 'DisplayName', 'LF');hold on; axis tight;
    plot(selectedFrameTimeStamps(timePeriod), y(timePeriod), '-r', 'DisplayName', 'RF');
    temp = aPaw{1} - aPaw{2};
    phaseDiff = abs(mod(temp + pi, 2*pi) - pi);
    plot(selectedFrameTimeStamps(timePeriod), phaseDiff(timePeriod), '-k', 'DisplayName', 'Diff', 'LineWidth', 2);
    legend({'LF', 'RF', 'Diff'})
    xlabel('time(second)');
    ylabel('phase (diff) after hilbert transformation');
    
    
    subplot(3, 1, 2)
    x = paw{1}; y = aPaw{1};
    plot(selectedFrameTimeStamps(timePeriod), x(timePeriod)/20 + 3, '-b'); axis tight; hold on
    plot(selectedFrameTimeStamps(timePeriod), y(timePeriod), '-r');
    legend({'LF raw', 'LF hilbert transform'});
    xlabel('time(second)');
    
    
    subplot(3, 1, 3)
    x = paw{2}; y = aPaw{2};
    plot(selectedFrameTimeStamps(timePeriod), x(timePeriod)/20 + 3, '-b'); axis tight; hold on
    plot(selectedFrameTimeStamps(timePeriod), y(timePeriod), '-r');
    legend({'LF raw', 'LF hilbert transform'});
    xlabel('time(second)');
    
    
end

%% body angles

session = '191007_003';
unit_id = 24;


disp('loading neuralData.mat...');
sessionFolder = fullfile(getenv('OBSDATADIR'), 'sessions', session);
load(fullfile(sessionFolder, 'neuralData.mat'));
spkTimes = spkTimes';
unit_spkRates = spkRates(find(unit_ids == unit_id), :);
cellNum = find(unit_ids == unit_id);


disp('loading runAnalyzed.mat...')
sessionFolder = fullfile(getenv('OBSDATADIR'), 'sessions', session);
load(fullfile(sessionFolder, 'runAnalyzed.mat'));



% plot zscores for body angles
ind(1) = knnsearch(frameTimeStamps, timeStamps(1));
ind(2) = knnsearch(frameTimeStamps, timeStamps(end));
frameTimeStamps = frameTimeStamps(ind(1) : ind(2));
bodyAngles = bodyAngles(ind(1) : ind(2));

rows = 3;
cols = 4;
plotInd = 1;
figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
subplot(rows, cols, plotInd);
plot(frameTimeStamps, zscore(bodyAngles), '-');
xlabel('time (second)');
ylabel('Z-scores of body angles');
box off


plotInd = plotInd + 1;
subplot(rows, cols, plotInd);

signifInds = find(abs(zscore(bodyAngles))>=5);
signifTimes = frameTimeStamps(signifInds);
timeWindow = [-0.2, 0.2];

plotPSTH2(session, cellNum, signifTimes, {'xLims', timeWindow});
xlabel('time (0 = signif body angle)');

%% corr coef b/w phase and neural firing rate, for individual steps

s.stepPercentiles = [40 60]; % only include steps with durations in between these percentile limits

stanceBins = stanceBins(startTimeInd:endTimeInd, :);

phase = struct();
neuralA = struct();
positions = struct();
zpositions = struct();
for i = 1:4
    
    % get swing start and end times
    swingStartInds = find(diff(~stanceBins(:, i))==1);
    swingStartTimes = selectedFrameTimeStamps(swingStartInds(1:end-1));
    stanceEndTimes = selectedFrameTimeStamps(swingStartInds(2:end)-1);
    sessionEvents = cat(2, swingStartTimes, stanceEndTimes);
    sessionEvents = sessionEvents(~isnan(sum(sessionEvents,2)),:); % remove nan entries
    
    % only take steps in middle of duration distribution
    durations = diff(sessionEvents,1,2);
    durationLimits = prctile(durations, s.stepPercentiles);
    inds = find(durations>durationLimits(1) & durations<durationLimits(2));
    sessionEvents = sessionEvents(durations>durationLimits(1) & durations<durationLimits(2), :);
    % times{i} = sessionEvents;
    
    stanceEndInds = swingStartInds(inds+1)-1;
    swingStartInds = swingStartInds(inds);
    
    
    for j = 1:length(sessionEvents)
        x = selectedFrameTimeStamps(swingStartInds(j):stanceEndInds(j));
        positionsTemp = pawsXLocations(swingStartInds(j):stanceEndInds(j), i);
        yPositioinsTemp = pawsYLocations(swingStartInds(j):stanceEndInds(j), i);
        zPositionsTemp = pawsZLocations(swingStartInds(j):stanceEndInds(j), i);
        positionsTemp = positionsTemp - mean(positionsTemp);
        yPositioinsTemp = yPositioinsTemp - mean(yPositioinsTemp); 
        zPositionsTemp = zPositionsTemp - mean(zPositionsTemp);
        positions.(s.pawNames{i})(j, :) = interp1(1:length(positionsTemp), positionsTemp, linspace(1, length(positionsTemp), 100));
        zpositions.(s.pawNames{i})(j, :) = interp1(1:length(zPositionsTemp), zPositionsTemp, linspace(1, length(zPositionsTemp), 100));
        ypositions.(s.pawNames{i})(j, :) = interp1(1:length(yPositioinsTemp), yPositioinsTemp, linspace(1, length(yPositioinsTemp), 100));
        % positions = lowpass(positions, 5, 250);
        hPositions = hilbert(positionsTemp);
        y = angle(hPositions) + pi;
        phase.(s.pawNames{i})(j, :) = interp1(x, y, linspace(x(1), x(end), 100));
        
        bins = neuralTimes > sessionEvents(j, 1) & neuralTimes < sessionEvents(j, 2);
        nx = neuralTimes(bins);
        ny = unit_spkRates(bins);
        neuralA.(s.pawNames{i})(j, :) = interp1(nx, ny, linspace(nx(1), nx(end), 100));
    end
        
end

selected = datasample(1:length([phase.RF]), 32);
plotInd = 1;
figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on

for i = 1:length(selected)
     if plotInd > 16
        figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
        plotInd = 1;
     end   
    
    subplot(4, 4, plotInd);
    plot([phase.RF(selected(i),:)] - mean([phase.RF(selected(i),:)]));
    hold on;
    plot(([neuralA.RF(selected(i),:)] - mean([neuralA.RF(selected(i),:)]))/30);
    axis tight;
    plot(([zpositions.RF(selected(i), :)] - mean([zpositions.RF(selected(i), :)]))*1000);
    plot(([positions.RF(selected(i), :)] - mean([positions.RF(selected(i), :)]))*200);
    plot(([ypositions.RF(selected(i), :)] - mean([ypositions.RF(selected(i), :)]))*200);
  
    xlabel('swing start -> stance end');
    legend('RF limb phase', 'neural FR', 'RF limb z position', 'RF limb x position', 'RF limb y position');
    
    plotInd = plotInd + 1;
end

figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
plot(mean(phase.RF - mean(phase.RF, 2), 1), 'LineWidth', 1);
plot(mean(neuralA.RF - mean(neuralA.RF, 2), 1)/30, 'LineWidth', 1);
plot(mean(zpositions.RF - mean(zpositions.RF, 2), 1)*1000, 'LineWidth', 1);
plot(mean(positions.RF - mean(positions.RF, 2), 1)*200, 'LineWidth', 1);
plot(mean(ypositions.RF - mean(ypositions.RF, 2), 1)*700, 'LineWidth', 1);
xlabel('swing start -> stance end');
legend('RF limb phase', 'neural FR', 'RF limb z position', 'RF limb x position', 'RF limb y position');

figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
c = [];
samples = datasample(1:length(neuralA.(s.pawNames{paw})), 20);
for i = 1:length(neuralA.(s.pawNames{paw}))
    
    x = zpositions.(s.pawNames{paw})(i, :);
    y = neuralA.(s.pawNames{paw})(i, :);
    [c(i, :), lags] = xcorr(x - mean(x), y - mean(y), 30, 'normalized');
    
    if ismember(i, samples)
        xrange = linspace(-30, 30, size(c, 2));
        plot(xrange, c(i, :), '-', 'LineWidth', 1, 'Color', [1 0.77, 0.16, 0.4]); hold on
    end
    
end

plot(xrange, mean(c, 1), '-', 'LineWidth', 4, 'Color', [0 0.4470 0.7410]); hold on;
minind = find(mean(c, 1) == min(mean(c, 1))); 
xlimit = get(gca,'xlim');
ylimit = get(gca,'ylim');
text(xrange(minind), min(mean(c, 1)), ['min corr = ' num2str(min(mean(c, 1)))], 'FontSize', 16);
text(xrange(minind), ylimit(2) - 0.1*(ylimit(2) - ylimit(1)), ['min.x = ' num2str(xrange(minind))], 'FontSize', 16);

line([xrange(minind), xrange(minind)], [ylimit(1), ylimit(2)], 'Color', [1 0.77 0.16], 'LineWidth', 4);

ylabel('corr coef');
xlabel('% of one step cycle')
title('cross correlation b/w RF paw z position and neural firing rates', 'FontSize', 18);

%% plot continuous test, chunk data version

% get velocity for each frame
if size(velTimeStamps, 1) == 1
    velTimeStamps = velTimeStamps';
end
frameVel = nan(length(frameTimeStamps), 1);

% find corresponding velocity for every frame times 
temp = frameTimeStamps(:, 1);
[~,ind] = histc(temp,[-inf;(velTimeStamps(1:end-1)+velTimeStamps(2:end))/2;inf]);
frameTimeStamps(:, 2) = velTimeStamps(ind(:, 1));
frameVel = vel(ind(:, 1))';

[runningBins, runningStartInds, runningEndInds] = isRunning(frameVel);

% processing paws x y z positions 
[locations, features] = fixTracking(locationsTable, frameTimeStamps, pixelsPerM);
botPawInds = find(contains(features, 'paw') & contains(features, '_bot'));
topPawInds = find(contains(features, 'paw') & contains(features, '_top'));

locationsPaws = nan(size(locations,1), 3, 4);  % frameNum X xyz X pawNum
locationsPaws(:,1:2,:) = locations(:,:,botPawInds);
locationsPaws(:,3,:) = locations(:,2,topPawInds);
locationsPawsPixels = locationsPaws;  % copy version in pixel coordinates before making subsequent transformations
locationsPaws(:,2,:) = locationsPaws(:,2,:) - nosePos(2); % subtract midline from all y values
locationsPaws(:,3,:) = (wheelCenter(2)-wheelRadius) - locationsPaws(:,3,:); % flip z and set s.t. top of wheel is zero

locationsPaws = locationsPaws/pixelsPerM;

pawsXLocations = squeeze(locationsPaws(startTimeInd:endTimeInd, 1, :));
pawsYLocations = squeeze(locationsPaws(startTimeInd:endTimeInd, 2, :));
pawsZLocations = squeeze(locationsPaws(startTimeInd:endTimeInd, 3, :));

% processing tail x y z positions

botTailBins = contains(features, 'tail') & contains(features, '_bot');
topTailBins = contains(features, 'tail') & contains(features, '_top');
locationsTail = nan(size(locations,1), 3, 2);  % frameNum X xyz X tailbase/mid
locationsTail(:,1:2,:) = locations(:,:,botTailBins);
locationsTail(:,3,:) = locations(:,2,topTailBins);
locationsTail(:,2,:) = locationsTail(:,2,:) - nosePos(2); % subtract midline from all y values
locationsTail(:,3,:) = (wheelCenter(2)-wheelRadius) - locationsTail(:,3,:); % flip z and set s.t. top of wheel is zero

locationsTail = locationsTail/pixelsPerM;

locationsTail = locationsTail(startTimeInd:endTimeInd, :, :);

% -----------------------------------------------------------------------%
% Figures paw position
figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
plotInd = 1;
for i = 1:4
    plotInd = plotContinuousData(pawsZLocations(:, i), unit_spkRates, plotInd, logical(runningBins), frameTimeStamps, neuralTimes, ...
                                 {'dataName', [s.pawNames{i} ' paw z positions']});
end


figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
plotInd = 1;
for i = 1:4
    plotInd = plotContinuousData(pawsXLocations(:, i), unit_spkRates, plotInd, logical(runningBins), frameTimeStamps, neuralTimes, ...
                                 {'dataName', [s.pawNames{i} ' paw x positions']});
end


figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
plotInd = 1;
for i = 1:4
    plotInd = plotContinuousData(pawsYLocations(:, i), unit_spkRates, plotInd, logical(runningBins), frameTimeStamps, neuralTimes, ...
                                 {'dataName', [s.pawNames{i} ' paw y positions']});
end


% ----------------------------------------------------------------------- %
% Figures velocity + body angles + tail angles 

figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
plotInd = 1;
% velocity
plotInd = plotContinuousData(frameVel, unit_spkRates, plotInd, logical(runningBins), frameTimeStamps, neuralTimes, ...
                                 {'dataName', 'velocity (m/s)'});
% bodyAngle
plotInd = plotContinuousData(bodyAngles, unit_spkRates, plotInd, logical(runningBins), frameTimeStamps, neuralTimes, ...
                                 {'dataName', 'bodyAngles'});                             
                             
% tailAngle
tailBase = squeeze(locationsTail(:, :, 1));
tailMid = squeeze(locationsTail(:, :, 2));

tailBase = fillmissing(tailBase, 'spline');
tailMid = fillmissing(tailMid, 'spline');

tailAngles = (tailBase(:, 2) - tailMid(:, 2)) ./ (tailBase(:, 1) - tailMid(:, 1));

plotInd = plotContinuousData(tailAngles, unit_spkRates, plotInd, logical(runningBins), frameTimeStamps, neuralTimes, ...
                                 {'dataName', 'tailAngles'});


% body z angles
noseZPositions = locations(startTimeInd:endTimeInd, :, 7);
noseZPositions(:, 2) = (wheelCenter(2)-wheelRadius) - noseZPositions(:, 2);
noseZPositions = noseZPositions / pixelsPerM;

tailBaseTop = locations(startTimeInd:endTimeInd, :, 5);
tailBaseTop(:, 2) = (wheelCenter(2)-wheelRadius) - tailBaseTop(:, 2);
tailBaseTop = tailBaseTop / pixelsPerM;
tailBaseTop = fillmissing(tailBaseTop, 'spline');

bodyZAngles = (noseZPositions(:, 2) - tailBaseTop(:, 2)) ./ (noseZPositions(:, 1) - tailBaseTop(:, 1));

plotInd = plotContinuousData(bodyZAngles, unit_spkRates, plotInd, logical(runningBins), frameTimeStamps, neuralTimes, ...
                                 {'dataName', 'bodyZAngles'});



% ----------------------------------------------------------------------- %
% Figure velocity in x y z axis

figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
plotInd = 1;
frameFreq = 250; % Hz
xVel = nan(length(frameTimeStamps)-1, 4);
for i = 1:4
    xVel(:, i) = abs(diff(pawsXLocations(:, i))) / (1/frameFreq);
    plotInd = plotContinuousData(xVel(:, i), unit_spkRates, plotInd, logical(runningBins), frameTimeStamps(2:end), neuralTimes, ...
                                 {'dataName', [s.pawNames{i} ' paw x velocity (m/s)']});
end


figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
plotInd = 1;
yVel = nan(length(frameTimeStamps)-1, 4);
for i = 1:4
    yVel(:, i) = abs(diff(pawsYLocations(:, i))) / (1/frameFreq);
    plotInd = plotContinuousData(yVel(:, i), unit_spkRates, plotInd, logical(runningBins), frameTimeStamps(2:end), neuralTimes, ...
                                 {'dataName', [s.pawNames{i} ' paw y velocity (m/s)']});
end



figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
plotInd = 1;
zVel = nan(length(frameTimeStamps)-1, 4);
for i = 1:4
    zVel(:, i) = abs(diff(pawsZLocations(:, i))) / (1/frameFreq);
    plotInd = plotContinuousData(zVel(:, i), unit_spkRates, plotInd, logical(runningBins), frameTimeStamps(2:end), neuralTimes, ...
                                 {'dataName', [s.pawNames{i} ' paw z velocity (m/s)']});
end




% ----------------------------------------------------------------------- %
% Figure acceleration in x y z axis

figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
plotInd = 1;
xAcc = nan(length(frameTimeStamps)-2, 4);

for i = 1:4
    xAcc(:, i) = abs(diff(xVel(:, i))) / (1/frameFreq);
    plotInd = plotContinuousData(xAcc(:, i), unit_spkRates, plotInd, logical(runningBins), frameTimeStamps(3:end), neuralTimes, ...
                                 {'dataName', [s.pawNames{i} ' paw x Acceleration (m/s^2)']});
end



figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
plotInd = 1;
yAcc = nan(length(frameTimeStamps)-2, 4);

for i = 1:4
    yAcc(:, i) = abs(diff(yVel(:, i))) / (1/frameFreq);
    plotInd = plotContinuousData(yAcc(:, i), unit_spkRates, plotInd, logical(runningBins), frameTimeStamps(3:end), neuralTimes, ...
                                 {'dataName', [s.pawNames{i} ' paw y Acceleration (m/s^2)']});
end


figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
plotInd = 1;
zAcc = nan(length(frameTimeStamps)-2, 4);

for i = 1:4
    zAcc(:, i) = abs(diff(zVel(:, i))) / (1/frameFreq);
    plotInd = plotContinuousData(zAcc(:, i), unit_spkRates, plotInd, logical(runningBins), frameTimeStamps(3:end), neuralTimes, ...
                                 {'dataName', [s.pawNames{i} ' paw z Acceleration (m/s^2)']});
end



%% 

sessions = {'200116_000' '200117_000' '200131_000'};

for i = 1:length(sessions)
    plotEphysLocomotorScatters3(sessions{i});
end

%%

c = struct();
    for paw = 1:4  
        for i = 1:length(neuralA.(s.pawNames{paw}))
            x = positions.(s.pawNames{paw})(i, :);
            y = ypositions.(s.pawNames{paw})(i, :);
            z = zpositions.(s.pawNames{paw})(i, :);
            neuA = neuralA.(s.pawNames{paw})(i, :);
            [c.(s.pawNames{paw}).x(i, :), lags] = xcorr(x - mean(x), neuA - mean(neuA), 30, 'normalized');
            [c.(s.pawNames{paw}).y(i, :), lags] = xcorr(y - mean(y), neuA - mean(neuA), 30, 'normalized');
            [c.(s.pawNames{paw}).z(i, :), lags] = xcorr(z - mean(z), neuA - mean(neuA), 30, 'normalized');        
        end
    end


maxValue = struct();
for i = 1:4
    maxValue.(s.pawNames{i}).x = max(abs(c.(s.pawNames{i}).x), [], 2);
    maxValue.(s.pawNames{i}).y = max(abs(c.(s.pawNames{i}).y), [], 2);
    maxValue.(s.pawNames{i}).z = max(abs(c.(s.pawNames{i}).z), [], 2);
    maxValue.(s.pawNames{i}).x(isnan(maxValue.(s.pawNames{i}).x)) = [];
    maxValue.(s.pawNames{i}).y(isnan(maxValue.(s.pawNames{i}).y)) = [];
    maxValue.(s.pawNames{i}).z(isnan(maxValue.(s.pawNames{i}).z)) = [];
end
    
rows = 2;
cols = 4;
plotInd = 0;

figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
for i = 1:4
    plotInd = plotInd + 1;
    subplot(rows, cols, plotInd);
    x = maxValue.(s.pawNames{i}).x;
    y = maxValue.(s.pawNames{i}).y;
    z = maxValue.(s.pawNames{i}).z;
    [h(i), p(i)] = ttest2(x, z)
    histogram(x, 'FaceColor','c'); hold on
    histogram(z, 'FaceColor','r');
    title(['\bf Data Distribution ', s.pawNames{i}, ' paw'],'FontSize',15)
    legend({'max x corr coef','max z corr coef'},'FontSize',11,'Location','NorthWest')
    ylabel('\bf Frequency','FontSize',12)
    
    plotInd = plotInd + 1;
    subplot(rows, cols, plotInd);
    boxplot([x,z],'Notch','on','Labels',{'max x corr coef','max z corr coef'})
    title('\bf Mean Difference b/w max x and z corr coef','Fontsize',15)
    
    yLimit=get(gca,'ylim');
    xLimit=get(gca,'xlim');
    text(xLimit(1) + 0.02*(xLimit(2) - xLimit(1)), yLimit(2)- 0.02*(yLimit(2) - yLimit(1)), ['h = ', num2str(h(i))], ...
        'FontWeight', 'bold');
    text(xLimit(1) + 0.02*(xLimit(2) - xLimit(1)), yLimit(2) - 0.06*(yLimit(2) - yLimit(1)), ['p = ', num2str(p(i))], ...
         'FontWeight', 'bold');
     
end

%%

figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
bar(ks);
legend('KS2', 'KS1');
ylabel('Unit Number')
set(gca, 'XTickLabel',str)










    
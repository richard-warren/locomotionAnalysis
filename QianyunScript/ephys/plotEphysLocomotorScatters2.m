function plotEphysLocomotorScatters2(session, unit_id, opts)

% settings
s.folder = fullfile(getenv('OBSDATADIR'), 'figures', 'ephys', 'modifiedSteps');  % folder in which the plots will be saved
% s.firingRatePencitiles = [30, 30]; % only include modified steps with the lowest 30% FR and highest 30% FR
s.xcorrFrequency = 1000;  % unit in Hz
s.markerSize = 500;
s.pawNames = {'LH', 'LF', 'RF', 'RH'};
s.pawColors = hsv(4);
s.rows = 4;
s.cols = 4;

% reassign settings contained in opts
if exist('opts', 'var'); for i = 1:2:length(opts); s.(opts{i}) = opts{i+1}; end; end


% loading files
disp('loading runAnalyzed.mat...')
sessionFolder = fullfile(getenv('OBSDATADIR'), 'sessions', session);
load(fullfile(sessionFolder, 'runAnalyzed.mat'));

disp('loading kinData.mat...');
if exist(fullfile(sessionFolder, 'kinData.mat'))
    load(fullfile(sessionFolder, 'kinData.mat'), 'kinData', 'stanceBins')
else
    [kinData, stanceBins] = getKinematicData(session);
end

disp('loading neuralData.mat...');
load(fullfile(sessionFolder, 'neuralData.mat'));
spkTimes = spkTimes';
unit_spkRates = spkRates(find(unit_ids == unit_id), :);
unit_spkTimes = spkTimes{find(unit_ids == unit_id), :};

% get rid of nans in neural firing rates
startInd = find(~isnan(unit_spkRates), 1, 'first');
endInd = find(~isnan(unit_spkRates), 1, 'last');
unit_spkRates = unit_spkRates(startInd:endInd);
neuralTimes = timeStamps(startInd:endInd);

% restricting wheel times and vel by start and end of neural times
vel = getVelocity(wheelPositions, .02, 1000); % should have the same length as wheelTimes
ind(1) = knnsearch(wheelTimes', neuralTimes(1));
ind(2) = knnsearch(wheelTimes', neuralTimes(end));
velTimeStamps = wheelTimes(ind(1) : ind(2));
vel = vel(ind(1) : ind(2));

% restricting frame times and tracking data by start and end of neural times

locationsTable = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw.csv')); % get raw tracking data
[locations, features] = fixTracking(locationsTable, frameTimeStamps, pixelsPerM);

startTimeInd = knnsearch(frameTimeStamps, neuralTimes(1));
endTimeInd = knnsearch(frameTimeStamps, neuralTimes(end));
frameTimeStamps = frameTimeStamps(startTimeInd : endTimeInd);
pawsXLocations = locations(startTimeInd:endTimeInd, 1, 9:12);
pawsYLocations = locations(startTimeInd:endTimeInd, 2, 9:12);
pawsZLocations = locations(startTimeInd:endTimeInd, 2, 1:4);
% --------------------------------------------------------------- %

% Figure 1
disp('plotting vel + bodyAngles + tailPosition + ...')
fprintf('%s: plotting unit %i', session, unit_id)
figure('name', sprintf('%s - unit %i', session, unit_id), ...
    'color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on


% velocity plot
plotInd = 1;
plotInd = plotContinuousData(vel, unit_spkRates, plotInd, {'dataName', 'velocity (m/s)'});


% bodyAngles
[bodyAnglesInterp, ~, unit_spkRatesInterp, ~] = dataInterpolation(bodyAngles, frameTimeStamps, unit_spkRates, neuralTimes);
plotInd = plotContinuousData(bodyAnglesInterp, unit_spkRatesInterp, plotInd, {'dataName', 'bodyAngles'});


% tail angles (xy plane)
tailBase_bot = locations(startTimeInd:endTimeInd, :, 13);
tailMid_bot = locations(startTimeInd:endTimeInd, :, 14);

inds1 = intersect(find(~isnan(tailBase_bot(:, 1))), find(~isnan(tailBase_bot(:, 2))));
inds2 = intersect(find(~isnan(tailMid_bot(:, 1))), find(~isnan(tailMid_bot(:, 2))));
inds = intersect(inds1, inds2);
tailBase_bot = tailBase_bot(inds, :);
tailMid_bot = tailMid_bot(inds, :);
varTimes = frameTimeStamps(inds, 1);

tailAngles = (tailBase_bot(:, 2) - tailMid_bot(:, 2))./(tailBase_bot(:, 1) - tailMid_bot(:, 1));

[tailAnglesInterp, ~, unit_spkRatesInterp, ~] = dataInterpolation(tailAngles, varTimes, unit_spkRates, neuralTimes);
plotInd = plotContinuousData(tailAnglesInterp, unit_spkRatesInterp, plotInd, {'dataName', 'tailAngles (xy plane)'});


% z position of the body
clear inds1 inds2 inds
nose_top = locations(startTimeInd:endTimeInd, :, 7);
tailBase_top = locations(startTimeInd:endTimeInd, :, 5);

inds1 = intersect(find(~isnan(nose_top(:, 1))), find(~isnan(nose_top(:, 2))));
inds2 = intersect(find(~isnan(tailBase_top(:, 1))), find(~isnan(tailBase_top(:, 2))));
inds = intersect(inds1, inds2);
nose_top = nose_top(inds, :);
tailBase_top = tailBase_top(inds, :);
varTimes = frameTimeStamps(inds, 1);

bodyZAngles = (nose_top(:, 2) - tailBase_top(:, 2))./(nose_top(:, 1) - tailBase_top(:, 1));

[bodyZAnglesInterp, ~, unit_spkRatesInterp, ~] = dataInterpolation(bodyZAngles, varTimes, unit_spkRates, neuralTimes);
plotInd = plotContinuousData(bodyZAnglesInterp, unit_spkRatesInterp, plotInd, {'dataName', 'bodyZAngles (xz plane)'});



% ----------------------------------------------------------------------- %
% Figure 2 paw positions

pawsXLocations = locations(startTimeInd:endTimeInd, 1, 9:12);
pawsYLocations = locations(startTimeInd:endTimeInd, 1, 1:4);
pawsZLocations = locations(startTimeInd:endTimeInd, 2, 1:4);


disp('plotting paw x position ...')
fprintf('%s: plotting unit %i', session, unit_id)
figure('name', sprintf('%s - unit %i', session, unit_id), ...
    'color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
plotInd = 1;

for i = 1:4
    temp = squeeze(pawsXLocations(:, 1, i));
    if any(isnan(temp))
        temp = fillmissing(temp, 'linear');
    end
    
    [tempInterp, ~, unit_spkRatesInterp, ~] = dataInterpolation(temp, frameTimeStamps, unit_spkRates, neuralTimes);
    plotInd = plotContinuousData(tempInterp, unit_spkRatesInterp, plotInd, {'dataName', [s.pawNames{i}, ' x position']});
end


disp('plotting paw y position...')
fprintf('%s: plotting unit %i', session, unit_id)
figure('name', sprintf('%s - unit %i', session, unit_id), ...
    'color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
plotInd = 1;

for i = 1:4
    temp = squeeze(pawsYLocations(:, 1, i));
    if any(isnan(temp))
        temp = fillmissing(temp, 'linear');
    end
    
    [tempInterp, ~, unit_spkRatesInterp, ~] = dataInterpolation(temp, frameTimeStamps, unit_spkRates, neuralTimes);
    plotInd = plotContinuousData(tempInterp, unit_spkRatesInterp, plotInd, {'dataName', [s.pawNames{i}, ' y position']});
end


disp('plotting paw z position...')
fprintf('%s: plotting unit %i', session, unit_id)
figure('name', sprintf('%s - unit %i', session, unit_id), ...
    'color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
plotInd = 1;

for i = 1:4
    temp = squeeze(pawsZLocations(:, 1, i));
    if any(isnan(temp))
        temp = fillmissing(temp, 'linear');
    end
    
    [tempInterp, ~, unit_spkRatesInterp, ~] = dataInterpolation(temp, frameTimeStamps, unit_spkRates, neuralTimes);
    plotInd = plotContinuousData(tempInterp, unit_spkRatesInterp, plotInd, {'dataName', [s.pawNames{i}, ' z position']});
    
end



disp('plotting paw x position ...')
fprintf('%s: plotting unit %i', session, unit_id)
figure('name', sprintf('%s - unit %i', session, unit_id), ...
    'color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
plotInd = 1;

for i = 1:4
    xVar = s.pawNames{i};
    paw.(xVar) = pawsXLocations(:, :, i) - mean(pawsXLocations(:, :, i));
    paw.(xVar) = lowpass(paw.(xVar), 5, 250);
    hPaw.(xVar) = hilbert(squeeze([paw.(xVar)]));
    aPaw.(xVar) = angle([hPaw.(xVar)]) + pi;
    
    [tempInterp, ~, unit_spkRatesInterp, ~] = dataInterpolation([aPaw.(xVar)], frameTimeStamps, unit_spkRates, neuralTimes);
    plotInd = plotContinuousData(tempInterp, unit_spkRatesInterp, plotInd, {'dataName', [xVar, ' phase']});    
end




%% interlimb coordination

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
pawName = {'LH', 'LF', 'RF', 'RH'};
figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize'));
plotInd = 0;
phaseDiff = [];
validStartInds = [];
validEndInds = [];

for i = 1:length(runningStartInds)
    paw = struct();
    hPaw = struct();
    aPaw = struct();
    trialPhaseDiff = [];
    startInd = [];
    endInd = [];
    ind = [];
    
    for pawInd = 1:4
        xVar = pawName{pawInd};
        paw.(xVar) = pawsXLocations(runningStartInds(i):runningEndInds(i), :, pawInd) - mean(pawsXLocations(runningStartInds(i):runningEndInds(i), :, pawInd));
        paw.(xVar) = lowpass(paw.(xVar), 5, 250);
        hPaw.(xVar) = hilbert(squeeze([paw.(xVar)]));
        aPaw.(xVar) = angle([hPaw.(xVar)]) + pi;
    end
    
    limb1 = 'LF';
    limb2 = 'RF';
    
    
    ind = vertcat(find(diff([aPaw.(limb1)]) < 0), find(diff([aPaw.(limb2)]) < 0)); 
    ind = sort(ind);
    startInd = ind(find(diff(ind)>10, 1, 'first') + 2);
    endInd = ind(find(diff(ind) > 10, 1, 'last'));
    
    temp = [aPaw.(limb1)] - [aPaw.(limb2)];

    chunkPhaseDiff = abs(mod(temp(startInd:endInd) + pi, 2*pi) - pi);
    phaseDiff = vertcat(phaseDiff, chunkPhaseDiff);
    validStartFrameInds(i, 1) = runningStartInds(i) + startInd - 1;
    validEndFrameInds(i, 1) = runningStartInds(i) + endInd - 1;
    
    if any(chunkPhaseDiff < 2)
        plotInd = plotInd + 1;
        
        subplot(4, 4, plotInd);
        plot([aPaw.(limb1)(startInd:endInd)], '-b', 'DisplayName', 'LF'); hold on
        plot([aPaw.(limb2)(startInd:endInd)] , '-r', 'DisplayName', 'RF');
        plot(chunkPhaseDiff, '-k', 'DisplayName', 'Diff', 'LineWidth', 2);
        legend({['phase of ', limb1], ['phase of ', limb2], 'phase Diff'})
        xlabel('frames');
        ylabel('phase (diff)');
        axis tight
        
        if plotInd >= 16
            figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize'));
            plotInd = 0;
        end
    end    
end














figure('color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize'));
plotInd = 0;
for i = 1:length(obsOnTimes)
    trialVel = vel(velTimeStamps >= obsOnTimes(i) & velTimeStamps <= obsOffTimes(i));
 
    if nanmean(trialVel) > 0.3
        
        ind(1) = find(frameTimeStamps >= obsOnTimes(i), 1, 'first');
        ind(2) = find(frameTimeStamps >= obsOffTimes(i), 1, 'first');
        for i = 1:4
            xVar = pawName{i};
            paw.(xVar) = pawsXLocations(ind(1):ind(2), :, i) - mean(pawsXLocations(ind(1):ind(2), :, i));
            paw.(xVar) = lowpass(paw.(xVar), 5, 250);
            hPaw.(xVar) = hilbert(squeeze([paw.(xVar)]));
            aPaw.(xVar) = angle([hPaw.(xVar)]) + pi;
        end
        
        limb1 = 'LF';
        limb2 = 'RF';
        
        temp = [aPaw.(limb1)] - [aPaw.(limb2)];
        trialPhaseDiff = abs(mod(temp + pi, 2*pi) - pi);
        
        if randn > 0.5
            plotInd = plotInd + 1;
            
            subplot(4, 4, plotInd);
            plot([aPaw.(limb1)], '-b', 'DisplayName', 'LF'); hold on
            plot([aPaw.(limb2)] , '-r', 'DisplayName', 'RF');
            plot(trialPhaseDiff, '-k', 'DisplayName', 'Diff', 'LineWidth', 2);
            legend({['phase of ', limb1], ['phase of ', limb2], 'phase Diff'})
            xlabel('frames');
            ylabel('phase (diff)');
            axis tight
            
            if plotInd >= 16
                break
            end
        end
        
    end
    
    phaseDiff = vertcat(phaseDiff, trialPhaseDiff);
    
end












end

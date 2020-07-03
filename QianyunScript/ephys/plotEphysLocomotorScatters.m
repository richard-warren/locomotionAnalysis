function plotEphysLocomotorScatters(session, unit_id, opts)

% settings
s.folder = fullfile(getenv('OBSDATADIR'), 'figures', 'ephys', 'modifiedSteps');  % folder in which the plots will be saved
% s.firingRatePencitiles = [30, 30]; % only include modified steps with the lowest 30% FR and highest 30% FR
s.xcorrFrequency = 1000;  % unit in Hz
s.markerSize = 500;
s.pawNames = {'LH', 'LF', 'RF', 'RH'};
s.pawColors = hsv(4);
s.rows = 3;
s.cols = 3;



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



% get experiment data
disp('getting experiment data...')
data = getExperimentData(session, {'isLeading','isTrialSuccess', 'isBigStep', 'obsHgt', 'firstModPaw'...
                                   'stepOverMaxHgt', 'controlStepHgt', 'noObsStepHgt', 'stepOverLength',  'controlStepLength'...
                                   'stepOverKinInterp', 'controlStepKinInterp'});
Data = data.data.sessions.trials;
flat = flattenData(Data, {'trial', 'paw', 'isLeading', 'isTrialSuccess', 'isBigStep', 'obsHgt', 'firstModPaw' ...
                          'stepOverMaxHgt', 'controlStepHgt', 'noObsStepHgt', 'stepOverLength',  'controlStepLength',...
                          'stepOverKinInterp', 'controlStepKinInterp'});
                      

% plot scatters of step heights and step length for 4 different limbs
% across 3 different step types

% getting step times
disp('getting step times...');
pawName = {'LH', 'LF', 'RF', 'RH'};
for paw = 1:4
    pawVar = pawName{paw};
    stepOverTimes.(pawVar) = getStepTimes(session, paw, 'stepOver', kinData);
    
    controlStepTimes.(pawVar) = getStepTimes(session, paw, 'control', kinData);
    
    noObsStepTimes.(pawVar) = getStepTimes(session, paw, 'noObs', kinData);
    
end


% getting step height for control steps and no obs steps
trialNum = length(kinData);
controlStepHgt = nan(1, 2);
noObsStepHgt = nan(1, 3);

for i = 1:trialNum
    for paw = 1:4
        if ~kinData(i).isTrialAnalyzed
            ind = find([flat.trial] == i & [flat.paw] == paw);
            flat(ind).rawControlStepHgt = nan;
            flat(ind).rawNoObsStepHgt = nan;
            flat(ind).controlStepPawInd = nan;
            flat(ind).noObsStepPawInd = nan;
        else   
            pawInd = [];
            for stepNum = 1:2                
                controlStepHgt(stepNum) = cell2mat(cellfun(@(x) max(squeeze(x(stepNum,3,:))),...
                                                   kinData(i).controlLocationsInterp(paw), 'UniformOutput', false));
                pawInd(stepNum) = paw;
            end
            ind = find([flat.trial] == i & [flat.paw] == paw);
            flat(ind).rawControlStepHgt = controlStepHgt;
            flat(ind).controlStepPawInd = pawInd;
            
            pawInd = [];
            for stepNum = 1:3
                noObsStepHgt(stepNum) = cell2mat(cellfun(@(x) max(squeeze(x(stepNum,3,:))),...
                                                   kinData(i).noObsLocationsInterp(paw), 'UniformOutput', false));
                pawInd(stepNum) = paw;
            end
            flat(ind).rawNoObsStepHgt = noObsStepHgt;
            flat(ind).noObsStepPawInd = pawInd;
        end
        
    end
end

% getting mean FR for different steps
disp('getting mean FR for steps...');

for i = 1:trialNum
    
    for paw = 1:4
        if ~kinData(i).isTrialAnalyzed
            ind = find([flat.trial] == i & [flat.paw] == paw);
            flat(ind).stepOverStepMeanFR = nan;
            flat(ind).controlStepMeanFR = nan;
            flat(ind).noObsStepMeanFR = nan;
        else                                    
            % meanFR for step over steps
            pawVar = pawName{paw};
            temp = unit_spkRates(timeStamps > stepOverTimes.(pawVar)(i).startTimes & timeStamps < stepOverTimes.(pawVar)(i).endTimes);
            ind = find([flat.trial] == i & [flat.paw] == paw);
            if any(isnan(temp))
                flat(ind).stepOverStepMeanFR = nan;
            else
                flat(ind).stepOverStepMeanFR = mean(temp);
            end
            
            % meanFR for control steps
            for j = 1:length(controlStepTimes.(pawVar)(i).startTimes)                
                temp = unit_spkRates(timeStamps > controlStepTimes.(pawVar)(i).startTimes(j) & timeStamps < controlStepTimes.(pawVar)(i).endTimes(j));
                ind = find([flat.trial] == i & [flat.paw] == paw);
                if any(isnan(temp))
                    flat(ind).controlStepMeanFR(j) = nan;
                else
                    flat(ind).controlStepMeanFR(j) = mean(temp);
                end
                
            end
            
            % meanFR for noOns steps
            for j = 1:length(noObsStepTimes.(pawVar)(i).startTimes)                
                temp = unit_spkRates(timeStamps > noObsStepTimes.(pawVar)(i).startTimes(j) & timeStamps < noObsStepTimes.(pawVar)(i).endTimes(j));
                ind = find([flat.trial] == i & [flat.paw] == paw);
                if any(isnan(temp))
                    flat(ind).noObsStepMeanFR(j) = nan;
                else
                    flat(ind).noObsStepMeanFR(j) = mean(temp);
                end                
            end
        end        
    end    
end
        
%%         
% plot scatter plots for steps!!
% Plot!!!
disp('plotting...')
fprintf('%s: plotting unit %i', session, unit_id)
figure('name', sprintf('%s - unit %i', session, unit_id), ...
    'color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on


xVar = {'stepOverMaxHgt', 'rawControlStepHgt', 'rawNoObsStepHgt'};
yVar = {'stepOverStepMeanFR', 'controlStepMeanFR', 'noObsStepMeanFR'};
conditionVar = {'paw', 'controlStepPawInd', 'noObsStepPawInd'};
plotInd = 1;
for i = 1:length(xVar)
    
    subplot(s.rows, s.cols, plotInd);
    [corrs, slopes] = scatterPlotRick([flat.(xVar{i})], [flat.(yVar{i})], [flat.(conditionVar{i})], ...
                                       {'scatAlpha', 0.2, 'conditionNames', {'LH', 'LF', 'RF', 'RH'}});
    xlabel(xVar{i});
    ylabel('mean FR');
    plotInd = plotInd + 1;
    
end


for paw = 1:4
    subplot(s.rows, s.cols, plotInd);
    bins = [flat.paw] == paw;
    temp = flat(bins);
    steps = cat(2, [temp.stepOverMaxHgt], [temp.rawControlStepHgt], [temp.rawNoObsStepHgt]);
    meanFR = cat(2, [temp.stepOverStepMeanFR], [temp.controlStepMeanFR], [temp.noObsStepMeanFR]);
    conditions = cat(2, repmat(1, 1, length([temp.paw])), repmat(2, 1, length([temp.controlStepPawInd])),...
                     repmat(3, 1, length([temp.noObsStepPawInd])));
    
    [corrs, slopes] = scatterPlotRick(steps, meanFR, conditions,...
                                      {'scatAlpha', 0.2, 'conditionNames', {'stepOverSteps', 'controlSteps', 'noObsSteps'}});
    hold on
    bins = ~isnan(steps) & ~isnan(meanFR);
    fit = polyfit(steps(bins), meanFR(bins), 1);
    plot(steps(bins), polyval(fit, steps(bins)), 'linewidth', 4, 'color', [0.8 0.6 0.93 1], 'DisplayName', 'All');
                                  
                                       
    xlabel([pawName{paw}, ' step Height']);
    ylabel('mean FR');
    plotInd = plotInd + 1;
end
    
 %%                     
% plot cros correlation b/w velocity and neural activity
unit_spkRates = unit_spkRates(~isnan(unit_spkRates));
timeStamps = timeStamps(~isnan(unit_spkRates));

vel = getVelocity(wheelPositions, .02, 1000); % should have the same length as wheelTimes
ind(1) = knnsearch(wheelTimes', timeStamps(1));
ind(2) = knnsearch(wheelTimes', timeStamps(end));
velTimeStamps = wheelTimes(ind(1) : ind(2));
vel = vel(ind(1) : ind(2));

% PLOT!
disp('plotting cross correlation...')
fprintf('%s: plotting unit %i', session, unit_id)
figure('name', sprintf('%s - unit %i', session, unit_id), ...
    'color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on
plotInd = 1;

if length(vel) == length(unit_spkRates)
    subplot(s.rows, s.cols, plotInd);
    [c, lags] = xcorr(vel', unit_spkRates', 10*s.xcorrFrequency);
    plot(lags, c, '-', 'LineWidth', 4);
    xlabel('time (seconds)');
    ylabel('crosscorr. b/w velocity and neural activities');
    box off
else
    display('length of velocity does not match length of unit firing rate!');
end


% plot zscores for body angles
ind(1) = knnsearch(frameTimeStamps, timeStamps(1));
ind(2) = knnsearch(frameTimeStamps, timeStamps(end));
frameTimeStamps = frameTimeStamps(ind(1) : ind(2));
bodyAngles = bodyAngles(ind(1) : ind(2));

plotInd = plotInd + 1;
subplot(s.rows, s.cols, plotInd);
plot(frameTimeStamps, zscore(bodyAngles), '-');
xlabel('time (second)');
ylabel('Z-scores of body angles');
box off



 % plot cros correlation b/w body angle and neural activity   
ind(1) = find(frameTimeStamps > timeStamps(1), 1, 'first');
ind(2) = find(frameTimeStamps < timeStamps(end), 1, 'last');
frameTimesInterp = frameTimeStamps(ind(1)): (1/s.xcorrFrequency) : frameTimeStamps(ind(2));
bodyAnglesInds = linspace(ind(1), ind(2), length(frameTimesInterp));
bodyAnglesInterp = interp1(1:length(bodyAngles), bodyAngles, bodyAnglesInds);

neuralTimeInds(1) = knnsearch(timeStamps', frameTimesInterp(1));
neuralTimeInds(2) = knnsearch(timeStamps', frameTimesInterp(end));
neuralTimes = timeStamps(neuralTimeInds(1):neuralTimeInds(2));
unit_spkRates = unit_spkRates(neuralTimeInds(1):neuralTimeInds(2));

if length(bodyAnglesInterp) == length(unit_spkRates)
    plotInd = plotInd + 1;
    subplot(s.rows, s.cols, plotInd);
    [c, lags] = xcorr(bodyAnglesInterp', unit_spkRates', 10*s.xcorrFrequency);
    plot(lags, c, '-', 'LineWidth', 4);
    xlabel('time (seconds)');
    ylabel('crosscorr. b/w body angles and neural activities');
    box off
else
    display('length of velocity does not match length of unit firing rate!');
end

%% interlimb coordination
disp('plotting interlimb correlation...')
fprintf('%s: plotting unit %i', session, unit_id)
figure('name', sprintf('%s - unit %i', session, unit_id), ...
    'color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on


% interlimb coordination
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'));
locationsTable = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw.csv')); % get raw tracking data
[locations, features] = fixTracking(locationsTable, frameTimeStamps, pixelsPerM);

pawsXLocations = locations(:, 1, 9:12);
pawName = {'LH', 'LF', 'RF', 'RH'};
paw = struct();
hPaw = struct();
aPaw = struct();


% data preparation

if size(wheelTimes, 1) == 1
    wheelTimes = wheelTimes';
end
frameVel = nan(length(frameTimeStamps), 1);
% find corresponding velocity for every frame times 
temp = frameTimeStamps(:, 1);
[~,ind] = histc(temp,[-inf;(wheelTimes(1:end-1)+wheelTimes(2:end))/2;inf]);
frameTimeStamps(:, 2) = wheelTimes(ind(:, 1));
frameVel = vel(ind(:, 1))';
% include only high speed times
timeInds =  find(frameVel > 0.3);

% perform the hilbert transform and get the phase estimation
for i = 1:4
    xVar = pawName{i};
    paw.(xVar) = pawsXLocations(timeInds, :, i) - mean(pawsXLocations(timeInds, :, i));
    paw.(xVar) = lowpass(paw.(xVar), 5, 250);
    hPaw.(xVar) = hilbert(squeeze([paw.(xVar)]));
    aPaw.(xVar) = angle([hPaw.(xVar)]) + pi;   
end

selectedTimes = frameTimeStamps(timeInds, 1);
exampleTimePeriod = 200000:201000;


% LF vs. RF
subplot(5, 4, 1:2);
limb1 = 'LF';
limb2 = 'RF';

temp = [aPaw.(limb1)] - [aPaw.(limb2)]; 
phaseDiff = abs(mod(temp + pi, 2*pi) - pi);

plot([aPaw.(limb1)(exampleTimePeriod)], '-b', 'DisplayName', 'LF');hold on; 
plot([aPaw.(limb2)(exampleTimePeriod)] , '-r', 'DisplayName', 'RF');
plot(phaseDiff(exampleTimePeriod), '-k', 'DisplayName', 'Diff', 'LineWidth', 2);
legend({['phase of ', limb1], ['phase of ', limb2], 'phase Diff'})
xlabel('frames');
ylabel('phase (diff)');
axis tight


subplot(5, 4, 3:4);
plot(selectedTimes, zscore(phaseDiff), 'o');
xlabel('time (second)');
ylabel('zscore of phase diff');


% LH vs. RH
subplot(5, 4, 5:6);
limb1 = 'LH'
limb2 = 'RH'

temp = [aPaw.(limb1)] - [aPaw.(limb2)]; 
phaseDiff = abs(mod(temp + pi, 2*pi) - pi);

plot([aPaw.(limb1)(exampleTimePeriod)], '-b', 'DisplayName', 'LF');hold on; 
plot([aPaw.(limb2)(exampleTimePeriod)] , '-r', 'DisplayName', 'RF');
plot(phaseDiff(exampleTimePeriod), '-k', 'DisplayName', 'Diff', 'LineWidth', 2);
legend({['phase of ', limb1], ['phase of ', limb2], 'phase Diff'})
xlabel('frames');
ylabel('phase (diff)');
axis tight



subplot(5, 4, 7:8);
plot(selectedTimes, zscore(phaseDiff), 'o');
xlabel('time (second)');
ylabel('zscore of phase diff');


% LF vs. LH
subplot(5, 4, 9:10);
limb1 = 'LF'
limb2 = 'LH'

temp = [aPaw.(limb1)] - [aPaw.(limb2)]; 
phaseDiff = abs(mod(temp + pi, 2*pi) - pi);

plot([aPaw.(limb1)(exampleTimePeriod)], '-b', 'DisplayName', 'LF');hold on; 
plot([aPaw.(limb2)(exampleTimePeriod)] , '-r', 'DisplayName', 'RF');
plot(phaseDiff(exampleTimePeriod), '-k', 'DisplayName', 'Diff', 'LineWidth', 2);
legend({['phase of ', limb1], ['phase of ', limb2], 'phase Diff'})
xlabel('frames');
ylabel('phase (diff)');
axis tight


subplot(5, 4, 11:12);
plot(selectedTimes, zscore(phaseDiff), 'o');
xlabel('time (second)');
ylabel('zscore of phase diff');


% LF vs. LH
subplot(5, 4, 13:14);
limb1 = 'RF'
limb2 = 'RH'

temp = [aPaw.(limb1)] - [aPaw.(limb2)]; 
phaseDiff = abs(mod(temp + pi, 2*pi) - pi);

plot([aPaw.(limb1)(exampleTimePeriod)], '-b', 'DisplayName', 'LF');hold on; 
plot([aPaw.(limb2)(exampleTimePeriod)] , '-r', 'DisplayName', 'RF');
plot(phaseDiff(exampleTimePeriod), '-k', 'DisplayName', 'Diff', 'LineWidth', 2);
legend({['phase of ', limb1], ['phase of ', limb2], 'phase Diff'})
xlabel('frames');
ylabel('phase (diff)');
axis tight


subplot(5, 4, 15:16);
plot(selectedTimes, zscore(phaseDiff), 'o');
xlabel('time (second)');
ylabel('zscore of phase diff');







end
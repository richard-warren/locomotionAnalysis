function plotEphysAndModifiedSteps(session, unit_id, paw, opts)

% settings
s.folder = fullfile(getenv('OBSDATADIR'), 'figures', 'ephys', 'modifiedSteps');  % folder in which the plots will be saved
s.firingRatePencitiles = [30, 30]; % only include modified steps with the lowest 30% FR and highest 30% FR
s.markerSize = 500;
s.pawNames = {'LH', 'LF', 'RF', 'RH'};
s.obsColor = [0.9 0.9 0.9; 0.3 0.74 0.93]; % colors for obstacle in the kinematics plot
s.colors = 'hsv'; % colors for conditions in the scatter plot

% reassign settings contained in opts
if exist('opts', 'var'); for i = 1:2:length(opts); s.(opts{i}) = opts{i+1}; end; end


% loading files
disp('loading runAnalyzed.mat...')
sessionFolder = fullfile(getenv('OBSDATADIR'), 'sessions', session);
load(fullfile(sessionFolder, 'runAnalyzed.mat'), 'frameTimeStamps', 'wheelPositions', 'wheelTimes' );

disp('loading kinData.mat...');
if exist(fullfile(sessionFolder, 'kinData.mat'))
    load(fullfile(sessionFolder, 'kinData.mat'), 'kinData', 'stanceBins')
else
    [kinData, stanceBins] = getKinematicData(session);
end

disp('loading neuralData.mat...');
load(fullfile(sessionFolder, 'neuralData.mat'));
spkTimes = spkTimes';
unit_spkTimes = spkTimes{find(unit_ids == unit_id), 1};  % an N-by-1 matrix
unit_spkRates = spkRates(find(unit_ids == unit_id), :);
clear spkRates spkTimes



% get experiment data
disp('getting experiment data...')
data = getExperimentData(session, {'isLeading','isTrialSuccess', 'isBigStep', 'obsHgt', 'firstModPaw'...
                                   'stepOverMaxHgt', 'controlStepHgt' 'isVentralContact', 'stepOverLength', 'stepOverKinInterp'});
Data = data.data.sessions.trials;
flat = flattenData(Data, {'trial', 'paw', 'isLeading', 'isTrialSuccess', 'isBigStep', 'obsHgt', 'firstModPaw' ...
                          'stepOverMaxHgt', 'controlStepHgt', 'isVentralContact', 'stepOverLength', 'stepOverKinInterp'});
                      
                  
flat = flat([flat.paw] == paw);


% calculate the mean firing rate during swing phase of step over steps for each trial and add it to the flat data
% structure
stepOverObsStepTimes = getStepOverObsStepTimes(session, paw);
controlStepTimes = getControlStepTimes(session, paw);


for i = 1:length(stepOverObsStepTimes)
     temp = unit_spkRates((timeStamps > stepOverObsStepTimes(i, 1)) & (timeStamps < stepOverObsStepTimes(i, 2)));
     if any(isnan(temp))
         flat(i).stepOverStepMeanFR = nan; 
     else
         flat(i).stepOverStepMeanFR = mean(temp);
     end        
end


for i = 1:length(controlStepTimes)
     temp = unit_spkRates((timeStamps > controlStepTimes(i, 1)) & (timeStamps < controlStepTimes(i, 2)));
     if any(isnan(temp))
         flat(i).controlStepMeanFR = nan; 
     else
         flat(i).controlStepMeanFR = mean(temp);
     end        
end
 

% calculate the mean velocity during swing phase of step over step for each
% trial and add it to the flate data structure

vel = getVelocity(wheelPositions, .02, 1000); % should have the same length as wheelTimes
velTimeStamps = wheelTimes;


for i = 1:length(stepOverObsStepTimes)
     temp = unit_spkRates((velTimeStamps > stepOverObsStepTimes(i, 1)) & (velTimeStamps < stepOverObsStepTimes(i, 2)));
     if any(isnan(temp))
         flat(i).meanVel = nan; 
     else
         flat(i).meanVel = mean(temp);
     end        
end
  

% calculate step max height and step length foor control steps
for i = 1:length(controlStepTimes)
    if kinData(i).isTrialAnalyzed
        temp = kinData(i).controlLocationsInterp(:, paw, :);
        temp = permute(cell2mat(temp), [2, 3, 1]);
        flat(i).controlStepMaxHgt = max(temp(3, :, 2));
        flat(i).controlStepLength = kinData(i).controlSwingLengths(end, paw);
    else
        flat(i).controlStepMaxHgt = nan;
        flat(i).controlStepLength = nan;
    end
end




% set binning variables and conditions 
bins = [flat.paw] == paw;
conditions = zeros(length(bins), 1);
conditions = conditions+1;

conditionNum = max(conditions);
if ischar(s.colors); s.colors = eval([s.colors '(conditionNum)']); end % set colorspace if color is specified as a string


% Plot!!!
fprintf('%s: plotting unit %i', session, unit_id)
figure('name', sprintf('%s - unit %i', session, unit_id), ...
            'color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on

subplot(1, 2, 1)
xVar = 'controlStepLength';        
[corrs, slopes] = scatterPlotRick([flat(bins).(xVar)], [flat(bins).meanFR], conditions(bins),...
                                   {'scatAlpha', 0.5});
xlabel([xVar, ' (meter)']);
ylabel('Mean FR During Swing Phase');

for i = 1:conditionNum
   text(max([flat(bins).(xVar)]), max([flat(bins).meanFR]) - (i-1)*5, ['condition', num2str(i), '  R^2 = ', num2str(corrs(i)^2)],...
        'Color', s.colors(i, :)); 
end
        
subplot(1, 2, 2)
xVar = 'controlStepMaxHgt';        
[corrs, slopes] = scatterPlotRick([flat(bins).(xVar)], [flat(bins).meanFR], conditions(bins),...
                                   {'scatAlpha', 0.5}); 
xlabel([xVar, ' (meter)']);
ylabel('Mean FR During Swing Phase');

for i = 1:conditionNum
   text(max([flat(bins).(xVar)]), max([flat(bins).meanFR]) - (i-1)*5,  ['condition', num2str(i), '  R^2 = ', num2str(corrs(i)^2)],...
        'Color', s.colors(i, :)); 
end
   




% plot the scatter histo plot

fprintf('%s: plotting unit %i', session, unit_id)
figure('name', sprintf('%s - unit %i', session, unit_id), ...
            'color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on


xVar = 'stepOverMaxHgt';        
scatterHistoRick([flat(bins).(xVar)], [flat(bins).stepOverStepMeanFR], 'groupId', conditions(bins), 'scatAlpha', 0.3, 'showCrossHairs', false,...
                 'groupHistoLineWidth', 3, 'xlabel', [xVar, ' (meter)'], 'ylabel', 'Mean FR During Swing Phase');


fprintf('%s: plotting unit %i', session, unit_id)
figure('name', sprintf('%s - unit %i', session, unit_id), ...
            'color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on            
xVar = 'stepOverLength'; 
scatterHistoRick([flat(bins).(xVar)], [flat(bins).meanFR], 'groupId', conditions(bins), 'scatAlpha', 0.3, 'showCrossHairs', false,...
                 'groupHistoLineWidth', 3, 'xlabel', [xVar, ' (meter)'], 'ylabel', 'Mean FR During Swing Phase');




             
% get the interpreted kinematics of the modified paw

bins = [flat.isLeading] == 1;

inds = find(bins == 1);
temp = 1;
stepOverKinInterp = nan(size([flat(1).stepOverKinInterp], 1), size([flat(1).stepOverKinInterp], 2), length(inds));
for i = inds
    stepOverKinInterp(:, :, temp) = [flat(i).stepOverKinInterp];
    temp = temp + 1;
end
stepOverKinInterp(2, :, :) = [];
trajectories = permute(stepOverKinInterp, [3, 1, 2]);


conditions = nan(length([flat(bins).meanFR]), 1);
FR_lowThresh = prctile([flat(bins).meanFR], 30);
FR_highThresh = prctile([flat(bins).meanFR], 70);
inds = find([flat(bins).meanFR] <= FR_lowThresh);
conditions(inds) = 1;
inds = find([flat(bins).meanFR] >= FR_highThresh);
conditions(inds) = 2;


conditionNum = max(conditions);
if ischar(s.colors); s.colors = eval([s.colors '(conditionNum)']); end % set colorspace if color is specified as a string


% PLOT!

fprintf('%s: plotting unit %i', session, unit_id)
figure('name', sprintf('%s - unit %i', session, unit_id), ...
            'color', 'white', 'units', 'pixels', 'position', get(0,'ScreenSize')); hold on

plotKinematics(trajectories, [flat(bins).obsHgt], conditions, ...
               'trialAlpha', .8, 'lineAlpha', 0, 'yLimZero', false, 'obsColors', s.obsColor)

t1 = text(0, -0.01, '\bf Binning variables: isLeading == 1');
t2 = text(0, -0.015, '\bf Steps with lowest 30% FR', 'Color', s.colors(1, :));
t3 = text(0, -0.02, '\bf Steps with highest 30% FR', 'Color', s.colors(2, :));






















end
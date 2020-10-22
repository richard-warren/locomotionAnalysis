%% compute experiment data from scratch

sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', 'whiskerTrimNotes');
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);  % remove not included sessions and empty lines
mice = unique(sessionInfo.mouse);

data = cell(1,length(mice));
parfor i=1:length(mice); data{i} = getExperimentData(sessionInfo(strcmp(sessionInfo.mouse, mice{i}),:), 'all'); end
data{1}.data = cellfun(@(x) x.data, data); data = data{1};
fprintf('saving...'); save(fullfile(getenv('OBSDATADIR'), 'matlabData', 'whiskerTrim_data.mat'), 'data'); disp('data saved!')


%% initializations

miceToExclude = {'den12'};
% miceToExclude = {};
global_config;
% fprintf('loading... '); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'whiskerTrim_data.mat'), 'data'); disp('whiskerTrim data loaded!')
fprintf('loading... '); load(fullfile('C:\Users\richa\Desktop\', 'matlabData', 'whiskerTrim_data.mat'), 'data'); disp('whiskerTrim data loaded!')  % for working from home
data.data = data.data(~ismember({data.data.mouse}, miceToExclude));

vars.condition = struct('name', 'condition', ...
    'levels', {{'bilateral full', 'unilateral full', 'unilateral int1', 'unilateral int2', 'unilateral int3', 'unilateral deltaOnly'}}, ...
    'levelNames', {{'biFull', 'unFull', 'un1', 'un2', 'un3', 'delta'}});
% vars.condition = struct('name', 'condition', ...
%     'levels', {{'bilateral full', 'unilateral int3', 'unilateral deltaOnly'}}, ...
%     'levelNames', {{'biFull', 'uniPartial', 'deltaOnly'}});
vars.isBigStep = struct('name', 'isBigStep', 'levels', [0 1], 'levelNames', {{'little step', 'big step'}});

colors = repmat([51 204 255]/255, length(vars.condition.levels), 1) .* fliplr(linspace(0,1,length(vars.condition.levels)))';

% set conditionals (for getDvMatrix)
conditionals.isFore = struct('name', 'isFore', 'condition', @(x) x==1);
conditionals.isLeading = struct('name', 'isLeading', 'condition', @(x) x==1);
conditionals.isLittleStep = struct('name', 'isBigStep', 'condition', @(x) x==0);

%% all bar plots
figure('position', [421 584 824 543], 'color', 'white', 'menubar', 'none');

% velocity
subplot(2,3,1)
dv = getDvMatrix(data, 'trialVel', vars.condition, {'mouse'});
barFancy(dv, 'ylabel', 'velocity (m/s)', 'levelNames', {}, 'colors', colors, barProperties{:}, ...
    'comparisons', [1 2; 1 3; 1 4; 1 5; 1 6], 'test', 'ttest')

% tail height
subplot(2,3,2)
dv = getDvMatrix(data, 'tailHgt', vars.condition, {'mouse'})*1000;
barFancy(dv, 'ylabel', 'tail height (mm)', 'levelNames', {}, 'colors', colors, barProperties{:}, ...
    'comparisons', [1 2; 1 3; 1 4; 1 5; 1 6], 'test', 'ttest')

% trial angle contra
subplot(2,3,3)
dv = getDvMatrix(data, 'trialAngleContra', vars.condition, {'mouse'});
barFancy(dv, 'ylabel', 'body angle', 'levelNames', {}, 'colors', colors, barProperties{:}, ...
    'comparisons', [1 2; 1 3; 1 4; 1 5; 1 6], 'test', 'ttest')

% success rate bars
subplot(2,3,4)
dv = getDvMatrix(data, 'isTrialSuccess', vars.condition, {'mouse'});
barFancy(dv, 'ylabel', 'success rate', 'levelNames', {vars.condition.levelNames}, 'colors', colors, barProperties{:}, ...
    'comparisons', [1 2; 1 3; 1 4; 1 5; 1 6], 'test', 'ttest')

% forepaw height
subplot(2,3,5)
dv = getDvMatrix(data, 'preObsHgt', vars.condition, {'mouse'}, [conditionals.isFore]);
barFancy(dv, 'ylabel', 'forepaw height (mm)', 'levelNames', {vars.condition.levelNames}, 'colors', colors, barProperties{:}, ...
    'comparisons', [1 2; 1 3; 1 4; 1 5; 1 6], 'test', 'ttest')

% correlation
subplot(2,3,6)
dv = getSlopeMatrix(data, {'obsHgt', 'preObsHgt'}, vars.condition, {'mouse'}, {'session'}, [conditionals.isLeading; conditionals.isFore], 'corr'); % perform regression for each session, then average slopes across sessions for each mouse
barFancy(dv, 'ylabel', 'paw-obstacle correlation', 'levelNames', {vars.condition.levelNames}, 'colors', colors, barProperties{:}, ...
    'comparisons', [1 2; 1 3; 1 4; 1 5; 1 6], 'test', 'ttest')

saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'whiskerTrimResults'), 'svg');


%% decision making initialization
flat = flattenData(data, ...
    [m.predictorsAll, {'mouse', 'isModPawLengthened', 'modPawDeltaLength', 'isBigStep', 'isLightOn', 'isModPawContra', ...
    'modPawOnlySwing', 'isTrialSuccess', 'condition', 'modPawPredictedDistanceToObs', 'modPawDistanceToObs', ...
    'modPawKinInterp', 'preModPawKinInterp', 'firstModPaw', 'preModPawDeltaLength', 'modSwingContacts'}]);

%% heatmaps
plotDecisionHeatmaps(flat, 'condition', 'condition', 'levels', vars.condition.levels, 'normalize', 'col', 'xLims', [-20 15], ...
    'modSwingContactsMax', 0, 'deltaMin', 0, 'successOnly', false, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'avgMice', false, 'plotMice', false, 'colors', colors, 'outcome', 'isModPawLengthened', 'plotProbs', false, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'whiskerTrimHeatmaps'));

%% trials scatters
plotDecisionTrials(flat, 'condition', 'condition', 'levels', vars.condition.levels, 'outcome', 'isBigStep', ...
    'modSwingContactsMax', false, 'deltaMin', 0, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'rowColors', colors, 'xLims', [-.11 .06], 'obsColor', obsColor, 'showHistos', true, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'whiskerTrimDecisionKin'));

%% model accuracies
[accuracies, ~, temp] = plotModelAccuracies(flat, m.predictors, 'isModPawLengthened', 'modelTransfers', [], ...
    'weightClasses', true, 'condition', 'condition', 'levels', vars.condition.levels, 'kFolds', 10, ...
    'modSwingContactsMax', m.modSwingContactsMax, 'deltaMin', m.deltaMin, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'colors', colors, 'barProps', [barProperties {'YLim', [.2 1], 'comparisons', [2 4; 2 6; 2 8; 2 10; 2 12], 'constantEdgeColor', [.15 .15 .15]}], ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'whiskerTrim_models'));

%% landing position entropy
plotEntropies(flat, 'condition', 'condition', 'levels', vars.condition.levels, 'poolMice', false, ...
    'modSwingContactsMax', 0, 'deltaMin', 0, 'successOnly', m.successOnly, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'colors', colors, 'barProps', [barProperties, {'comparisons', [1 2; 1 3; 1 4; 1 5; 1 6]}], ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'whiskerTrim_entropy'));
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'whiskerTrim_histos'), 'svg')

%% landing distance
figure('position', [200 472.00 300 255], 'color', 'white', 'menubar', 'none');
dv = getDvMatrix(data, 'modPawDistanceToObsAbs', [vars.condition], {'mouse'}, [figConditionals; conditionals.modPawOnlySwing; conditionals.isLittleStep])*1000;
barFancy(dv, 'ylabel', 'little step landing distance (mm)', 'levelNames', {vars.condition.levelNames}, 'colors', repmat(colors,2,1), barProperties{:}, ...
    'comparisons', [1 2; 1 3; 1 4; 1 5; 1 6], 'test', 'ttest')
saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'whiskerTrim_landingDistance'), 'svg');

%% (temp; to compare ispi and contra paw first) heatmaps
plotDecisionHeatmaps(flat(strcmp({flat.condition}, 'unilateral full')), 'condition', 'isModPawContra', 'levels', {0,1}, 'normalize', 'col', 'xLims', [-20 15], ...
    'modSwingContactsMax', 0, 'deltaMin', 0, 'successOnly', false, 'modPawOnlySwing', m.modPawOnlySwing, 'lightOffOnly', m.lightOffOnly, ...
    'avgMice', true, 'plotMice', true, 'colors', colors, 'outcome', 'isModPawLengthened', 'plotProbs', false);






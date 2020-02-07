%% compute experiment data from scratch

sessionInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'experimentMetadata.xlsx'), 'Sheet', 'whiskerTrimNotes');
sessionInfo = sessionInfo(sessionInfo.include==1 & ~cellfun(@isempty, sessionInfo.session),:);  % remove not included sessions and empty lines
mice = unique(sessionInfo.mouse);

data = cell(1,length(mice));
parfor i=1:length(mice); data{i} = getExperimentData(sessionInfo(strcmp(sessionInfo.mouse, mice{i}),:), 'all'); end
data{1}.data = cellfun(@(x) x.data, data); data = data{1};
fprintf('saving...'); save(fullfile(getenv('OBSDATADIR'), 'matlabData', 'whiskerTrim_data.mat'), 'data'); disp('data saved!')


%% initializations

global_config;
fprintf('loading... '); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'whiskerTrim_data.mat'), 'data'); disp('whiskerTrim data loaded!')
% data.data = data.data(~ismember({data.data.mouse}, miceToExclude));

vars.condition = struct('name', 'condition', ...
    'levels', {{'bilateral full', 'unilateral full', 'unilateral int1', 'unilateral int2', 'unilateral int3', 'unilateral deltaOnly'}}, ...
    'levelNames', {{'biFull', 'unFull', 'un1', 'un2', 'un3', 'delta'}});
% vars.condition = struct('name', 'condition', ...
%     'levels', {{'bilateral full', 'unilateral int3', 'unilateral deltaOnly'}}, ...
%     'levelNames', {{'biFull', 'uniPartial', 'deltaOnly'}});

colors = repmat([51 204 255]/255, length(vars.condition.levels), 1) .* fliplr(linspace(0,1,length(vars.condition.levels)))';

%% decision making

% settings
modPawOnlySwing = true;  % if true, only include trials where the modified paw is the only one in swing
lightOffOnly = false;  % whether to restrict to light on trials
successOnly = false;  % whether to only include successful trials


% heatmaps
plotDecisionHeatmaps(data, 'condition', 'condition', 'levels', vars.condition.levels, ...
    'successOnly', successOnly, 'modPawOnlySwing', modPawOnlySwing, 'lightOffOnly', lightOffOnly, ...
    'avgMice', true, 'plotMice', true, 'colors', colors, 'binNum', 50, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'whiskerTrim_heatmaps'));

% trials scatters
plotDecisionTrials(data, 'condition', 'condition', 'levels', vars.condition.levels, ...
    'successOnly', successOnly, 'modPawOnlySwing', modPawOnlySwing, 'lightOffOnly', lightOffOnly, ...
    'colors', decisionColors, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'whiskerTrim_decisionKin'));

% model accuracies
plotModelAccuracies(data, 'condition', 'condition', 'levels', vars.condition.levels, ...
    'successOnly', successOnly, 'modPawOnlySwing', modPawOnlySwing, 'lightOffOnly', lightOffOnly, ...
    'colors', [colors; .2 .2 .2], 'barProps', barProperties, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'whiskerTrim_models'));

% decision threshold
plotDecisionThresholds(data, 'condition', 'condition', 'levels', vars.condition.levels, ...
    'successOnly', successOnly, 'modPawOnlySwing', modPawOnlySwing, 'lightOffOnly', lightOffOnly, ...
    'colors', colors, 'barProps', barProperties, ...
    'saveLocation', fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'whiskerTrim_thresholds'));








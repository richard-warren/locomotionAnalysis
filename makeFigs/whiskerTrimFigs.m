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
% fprintf('loading... '); load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'whiskerTrim_data.mat'), 'data'); disp('whiskerTrim data loaded!')
fprintf('loading... '); load(fullfile('C:\Users\richa\Desktop\', 'matlabData', 'whiskerTrim_data.mat'), 'data'); disp('whiskerTrim data loaded!')  % for working from home

vars.condition = struct('name', 'condition', ...
    'levels', {{'bilateral full', 'unilateral full', 'unilateral int1', 'unilateral int2', 'unilateral int3', 'unilateral deltaOnly'}}, ...
    'levelNames', {{'biFull', 'unFull', 'un1', 'un2', 'un3', 'delta'}});
% vars.condition = struct('name', 'condition', ...
%     'levels', {{'bilateral full', 'unilateral int3', 'unilateral deltaOnly'}}, ...
%     'levelNames', {{'biFull', 'uniPartial', 'deltaOnly'}});

colors = repmat([51 204 255]/255, length(vars.condition.levels), 1) .* fliplr(linspace(0,1,length(vars.condition.levels)))';

% set conditionals (for getDvMatrix)
conditionals.isFore = struct('name', 'isFore', 'condition', @(x) x==1);
conditionals.isLeading = struct('name', 'isLeading', 'condition', @(x) x==1);

%% plot everything
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
barFancy(dv, 'ylabel', 'paw:obstacle correlation', 'levelNames', {vars.condition.levelNames}, 'colors', colors, barProperties{:}, ...
    'comparisons', [1 2; 1 3; 1 4; 1 5; 1 6], 'test', 'ttest')

saveas(gcf, fullfile(getenv('OBSDATADIR'), 'papers', 'hurdles_paper1', 'figures', 'matlabFigs', 'whiskerTrimResults'), 'svg');






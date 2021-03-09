function [data, cellInfo] = getClusteringData()

% create struct of tables containing per cell info on group importance,
% responses, and locations // use to explore functional clustering

% todo: refactor to not require aggregateResponses?


% inits
fprintf('getting clustering data... ')
load('E:\lab_files\paper2\modelling\response_aggregates.mat', 'aggregates', 'cellInfo')
groupInfo = readtable('C:\Users\richa\Desktop\github\locomotionAnalysis\paper2\glm\residual_predictorSettings.xlsx', ...
    'sheet', 'groups', 'ReadRowNames', true);
cellInfoCurrent = getUnitInfo();

% check that cellInfo is up to date!
if ~isequal(cellInfo(:, {'session', 'unit'}), cellInfoCurrent(:, {'session', 'unit'}))
    disp('WARNING! cellInfo is not up to date. Rerun aggregate responses please :)')
end

% load sample model to find group names
load(fullfile('E:\lab_files\paper2\modelling\glms\residual_glms\', ...
    sprintf('%s_cell_%i_glm.mat', cellInfo.session{1}, cellInfo.unit(1))), 'models');
groups = models.Properties.RowNames(2:end);

% make data storage struct (one entry per predictor group)
% each group gets its own table of responses and importances for each unit
xSmps = size(aggregates{1,'aggregate'}{1},2);  % number of x axis samples
data = struct();
for i = 1:length(groups)
    data.(groups{i}).tbl = table(cellInfo.session, cellInfo.unit, nan(height(cellInfo), xSmps), nan(height(cellInfo),1), ...
        'VariableNames', {'session', 'unit', 'response', 'importance'});
end

% get responses for each group
for i = 1:length(groups)
    predictor = groupInfo{groups{i}, 'response'}{1};  % representative predictor for group
    data.(groups{i}).predictor = predictor;
    xlims = aggregates{predictor, 'xLims'};
    data.(groups{i}).x = linspace(xlims(1), xlims(2), xSmps);
    data.(groups{i}).tbl.response = aggregates{predictor, 'aggregate'}{1};
end


% add columns to cellInfo for importance of each feature
inits = cell(1, length(groups));
inits(:) = {nan(height(cellInfo), 1)};
tbl = table(inits{:}, 'VariableNames', cellfun(@(x) ['importance_' x], groups, 'UniformOutput', false));
cellInfo = cat(2, cellInfo, tbl);

for i = 1:height(cellInfo)

    % load importances
    file = fullfile('E:\lab_files\paper2\modelling\glms\residual_glms\', ...
        sprintf('%s_cell_%i_glm.mat', cellInfo.session{i}, cellInfo.unit(i)));
    
    if exist(file, 'file')
        load(file, 'models');  % load sample models
        for j = 1:length(groups)
            if ismember(groups{j}, models.Properties.RowNames)
                data.(groups{j}).tbl.importance(i) = models{groups{j}, 'dev'};
                cellInfo{i, ['importance_' groups{j}]} = models{groups{j}, 'dev'};
            end
        end
    end
end


% save
save('E:\lab_files\paper2\modelling\clustering_data.mat', 'data', 'cellInfo')
fprintf('all done!\n')




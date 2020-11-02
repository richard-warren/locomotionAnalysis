function [data, cellInfo] = getClusteringData()

% create struct of tables containing per cell info on group importance,
% responses, and locations // use to explore functional clustering

% function plotGroupResponses()

% plot heatmaps of representative response variable for each response group
% // sort by deviance explained to sanity check deviance metric


% inits
fprintf('getting clustering data... ')
load('E:\lab_files\paper2\modelling\response_aggregates.mat', 'aggregates', 'cellInfo')
load(fullfile('E:\lab_files\paper2\modelling\glms\residual_glms\', ...
    [cellInfo.session{1} '_cell_' num2str(cellInfo.unit(1)) '_glm.mat']), 'models');  % load sample models
groupInfo = readtable('C:\Users\richa\Desktop\github\locomotionAnalysis\paper2\glm\predictorSettings.xlsx', ...
    'sheet', 'groups', 'ReadRowNames', true);
groups = models.Properties.RowNames(2:end);

% load ccf locations for all cells
files = dir('E:\lab_files\paper2\histo\registration\*_registration.mat');
registration = cell(1, length(files));
for i = 1:length(files)
    reg = load(fullfile(files(i).folder, files(i).name), 'registration');
    registration{i} = reg.registration;
end
registration = cat(1, registration{:});

% make data storage struct
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

% get importance for each unit
for i = 1:height(cellInfo)

    % load importances
    file = fullfile('E:\lab_files\paper2\modelling\glms\residual_glms\', ...
        sprintf('%s_cell_%i_glm.mat', cellInfo.session{i}, cellInfo.unit(i)));
    
    if exist(file, 'file')
        load(file, 'models');  % load sample models
        for j = 1:length(groups)
            if ismember(groups{j}, models.Properties.RowNames)
                data.(groups{j}).tbl.importance(i) = max(models{groups{j}, 'dev'}, 0);
            end
        end
    end
end

% get target locations
ephysInfo = readtable('Y:\loco\obstacleData\spreadSheets\ephysInfo.xlsx', ...
    'sheet', 'ephysInfo', 'ReadRowNames', true);
target = cell(height(cellInfo), 1);
for i = 1:length(target); target{i} = ephysInfo{cellInfo.session{i}, 'target'}{1}; end
cellInfo = cat(2, cellInfo, table(target, 'VariableNames', {'target'}));

% todo: get ccf locations
ccf = nan(height(cellInfo), 3);
for i = 1:height(cellInfo)
    bin = strcmp(registration.session, cellInfo.session{i}) & registration.unit==cellInfo.unit(i);
    if any(bin); ccf(bin,:) =  registration{bin, 'ccfMm'}; end
end
cellInfo = cat(2, cellInfo, table(ccf, 'VariableNames', {'ccf'}));

% save
save('E:\lab_files\paper2\modelling\clustering_data.mat', 'data', 'cellInfo')
fprintf('all done!\n')




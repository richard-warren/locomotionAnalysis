% function plotGroupResponses()

% plot heatmaps of representative response variable for each response group
% // sort by deviance explained to sanity check deviance metric


% inits
load('E:\lab_files\paper2\modelling\response_aggregates.mat', 'aggregates', 'cellInfo')
load(fullfile('E:\lab_files\paper2\modelling\glms\residual_glms\', ...
    [cellInfo.session{1} '_cell_' num2str(cellInfo.unit(1)) '_glm.mat']), 'models');  % load sample models
groupInfo = readtable('C:\Users\richa\Desktop\github\locomotionAnalysis\paper2\glm\predictorSettings.xlsx', ...
    'sheet', 'groups', 'ReadRowNames', true);
groups = models.Properties.RowNames(2:end);
responses = groupInfo{groups, 'response'};  % predictors to associate with each group

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

%% get importance for each unit
for i = 1:height(cellInfo)
    disp(i)
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

%%

close all; figure('color', 'white', 'position', [852.00 2.00 428.00 1408.00])


% inits
traceNum = 10;
colors = lines(length(groups));

% group inits
i=7;
tbl = data.(groups{i}).tbl;
x = data.(groups{i}).x;
predictor = data.(groups{i}).predictor;
[~, sortInds] = sort(tbl.importance);


% histogram
subplot(5,1,1)
histogram(tbl.importance, 50, 'EdgeColor', 'none', 'FaceColor', colors(i,:))
xlabel('deviance explained')
set(gca, 'box', 'off')
title(groups{i})

% responses
subplot(5,1,2); hold on
plot(x, tbl.response(sortInds(1:traceNum),:)', 'color', [0 0 0 .4])
plot(x, tbl.response(sortInds(end-traceNum:end),:)', 'color', [colors(i,:) .6], 'LineWidth', 2)
ylabel('firing rate (z score)')
xlabel(predictor, 'Interpreter', 'none')
set(gca, 'box', 'off', 'XLim', [x(1) x(end)])

% heatmap
subplot(5,1,3:5)
imagesc(x, 1:height(tbl), tbl.response(sortInds, :));
colormap gray
xlabel(predictor, 'Interpreter', 'none')
set(gca, 'box', 'off', 'ytick', [], 'TickDir', 'out', 'ydir', 'normal')









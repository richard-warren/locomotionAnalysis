%% general 'scratch pad' for modelling paper

sessions = getEphysSessions();

%% perform various analyses on all sessions

sessions = getEphysSessions;
% sessions = sessions(1:33);  % temp

overwrite = false;
tic

% parpool('local', 4);  % set number of workers
for i = 1:length(sessions)
    folderSes = fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i});
    folder = fullfile(getenv('SSD'), 'modelling');
    
    try
        % format ephys data
        if overwrite || ~exist(fullfile(folderSes, 'neuralData.mat'), 'file')
            formatEphysData(sessions{i})
        end
        
        % predictors
        if overwrite || ~exist(fullfile(folder, 'predictors', [sessions{i} '_predictors.mat']), 'file')
            getPredictors(sessions{i}, 'plotPredictors', true, 'visible', 'off')
        end
        
        % neural responses
        if overwrite || ~exist(fullfile(folder, 'responses', [sessions{i} '_responses.mat']), 'file')
            getNeuralResponses(sessions{i})
        end
        
        % feature importance
        if overwrite || ~exist(fullfile(folder, 'importance', [sessions{i} '_importance.mat']), 'file')
            getFeatureImportance(sessions{i})
        end
        
        % plot neural responses
%         plotNeuralResponses(sessions{i}, 'visible', false, 'showImportance', true)
    
    catch exception
        fprintf('%s: PROBLEM! -> %s\n', sessions{i}, exception.identifier)
    end
end
toc


%% look into MI

% load(fullfile(getenv('SSD'), 'modelling', 'aggregates', 'aggregates.mat'), 'aggregates', 'cellInfo');
mi = cat(2, aggregates.mi{:})';  % (predictor X cells) matrix of mutual information

[X, removed] = rmmissing(mi);
names = aggregates.Properties.RowNames(~removed);
X = (X-mean(X,2));
X = X ./ vecnorm(X')';
corrMat = X * X';

%% fit gmms

% corrMat = corrMat(randperm(size(corrMat,1)),:);
corrMat = rmmissing(mi)';
corrMat = zscore(corrMat,0,2);

% pca
pcs = 10;
[coeff, score, ~, ~, explained] = pca(corrMat);
projection = coeff(:,1:pcs) * score(:,1:pcs)';  % projection onto first pcs
projection = projection' + mean(corrMat,1);

nGroups = 8;
gmm = fitgmdist(score(:,1:pcs), nGroups, ...
    'CovarianceType', 'full', 'SharedCovariance', false, ...
    'RegularizationValue', .01, 'Replicates', 10, 'Options', statset('MaxIter', 1000));
[groups, ~, ~, ~, maha] = gmm.cluster(score(:,1:pcs));
[~, sortInds] = sort(groups);


% groups = clusterdata(corrMat, 'maxclust', nGroups);
% [~, sortInds] = sort(groups);


colors = lines(nGroups);
close all; figure('Units', 'normalized', 'position', [0 0 1 1])
% subplot(1,3,1); imagesc(corrMat)  % original
% set(gca, 'YTick', 1:size(corrMat,1), 'yticklabel', names, ...
%     'xtick', 1:size(corrMat,2), 'XTickLabel', names, 'XTickLabelRotation', 90)

subplot(1,2,1); imagesc(corrMat(sortInds,:)); hold on  % original
scatter(zeros(1,size(corrMat,1)), 1:size(corrMat,1), 50, colors(groups(sortInds),:), 'filled')
set(gca, 'xlim', [0 size(corrMat,2)], 'xtick', 1:size(corrMat,2), 'XTickLabel', names, 'XTickLabelRotation', 90, 'TickDir', 'out')

subplot(1,2,2); imagesc(projection(sortInds,:)); hold on
scatter(zeros(1,size(corrMat,1)), 1:size(corrMat,1), 50, colors(groups(sortInds),:), 'filled')
set(gca, 'xlim', [0 size(corrMat,2)], 'xtick', 1:size(corrMat,2), 'XTickLabel', names, 'XTickLabelRotation', 90, 'TickDir', 'out')
% set(gca, 'xtick', 1:size(corrMat,2), 'XTickLabel', names, 'XTickLabelRotation', 90, 'TickDir', 'out')












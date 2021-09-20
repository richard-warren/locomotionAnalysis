%% compare models with different numbers of predictors, starting with null
% model (time basis only), then 'base' model containing things like wheel
% vel and lick times, then adding locomotion phase, and finally add all the
% fancy kinematic vars... idea is to test whether all the fancy stuff
% actually helps us predict anything



%% compute design matrices

data = getUnitInfo();
sessions = unique(data.session);
spreadsheet = fullfile(getenv('GITDIR'), 'locomotionAnalysis', 'paper2', 'glm', 'settings', 'sequential_predictorSettings.xlsx');

tic
parfor i = 1:length(sessions)
    try
        filename = fullfile(getenv('SSD'), 'paper2', 'modelling', ...
            'designMatrices', 'sequential', [sessions{i} '_designMatrix.mat']);
        makeDesignMatrix(sessions{i}, spreadsheet, ...
            'timeDegrees', 3, 'outputFileName', filename);
    catch exception
        fprintf('%s: PROBLEM! -> %s\n', sessions{i}, exception.identifier)
    end
end
fprintf('finished in %.1f minutes\n', toc/60)


%% fit model for single unit

% bad
% session = '191009_003'; unit = 267;
% session = '200713_000'; unit = 285;
session = '201012_000'; unit = 93;

% good
% session = '200116_000'; unit = 128;
% session = '200708_000'; unit = 118;

[models, fitdata] = fitSequentialGlms(session, unit);
% plotSequentialGlms(session, unit);


% plot all predictions
% load(['E:\lab_files\paper2\modelling\glms\sequential_glms\' session '_cell_' num2str(unit) '_glm.mat'], ...
%     'models', 'fitdata');
close all;
figure('color', 'white', 'position', [2.00 722.00 1278.00 634.00]); hold on;
plot(fitdata.t, [fitdata.yRate; fitdata.yhat'])
plot(fitdata.t, models{end, 'model'}{1}.fit_preval)

%% explore deviances
fullDeviances = model.glmnet_fit.dev;
partialModels = models{end, 'model'}{1}.models;
partialDeviances = nan(length(fullDeviances), length(partialModels));  % lambda X model
for i = 1:length(partialModels)
    partialDeviances(:, i) = partialModels{i}.dev;
end

close all;
figure('color', 'white', 'menubar', 'none', 'position', [454.00 841.00 560.00 420.00]); hold on
plot(partialDeviances, 'color', [0 0 0 .4])
plot(fullDeviances, 'color', lines(1), 'LineWidth', 2)
xlabel('\lambda')
plot(mean(partialDeviances, 2), 'color', 'red')
lambda_id = models{end, 'model'}{1}.lambda_min_id;
plot([lambda_id lambda_id], ylim)

%% fit models for all units

% settings
overwrite = true;


data = getUnitInfo();
tic; fprintf('\nfitting sequential GLMs for %i neurons...\n', height(data))

parfor i = 1:height(data)  % individual sessions are repeated for each neuron in session
    session = data{i, 'session'}{1};
    unit = data{i, 'unit'};
    
    try
        % fit models
        filename = fullfile(getenv('SSD'), 'paper2', 'modelling', 'glms', 'sequential_glms', ...
            [session '_cell_' num2str(unit) '_glm.mat']);
        if overwrite || ~exist(filename, 'file')
            fprintf('(%3i/%i) %s, cell %3i: fitting GLMs\n', i, height(data), session, unit);
            fitSequentialGlms(session, unit, 'verbose', false, 'parallel', false, 'save', true);
        end

        % plot (temp commented out)
        plotFilename = fullfile(getenv('SSD'), 'paper2', 'modelling', 'glms', 'sequential_glms', ...
            [session '_cell_' num2str(unit) '_glm.png']);
        if (overwrite || ~exist(plotFilename, 'file')) && exist(filename, 'file')
            plotSequentialGlms(session, unit, 'save', true, 'visible', false)
        end
    catch exception
        fprintf('%s: PROBLEM! -> %s, unit %i\n', session, unit, exception.identifier)
    end
end
fprintf('\nfinished in %.1f minutes\n', toc/60)

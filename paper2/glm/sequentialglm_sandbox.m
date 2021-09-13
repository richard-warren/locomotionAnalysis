%% description goes here!


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

fitSequentialGlms('201227_000', 194)


%%  STUFF IS NOT FIXED BEYOND THIS POINT!

%% fit models for all units

% settings
overwrite = true;


data = getUnitInfo();
tic; fprintf('\nfitting residual GLMs for %i neurons...\n', height(data))

parfor i = 1:height(data)  % individual sessions are repeated for each neuron in session
    session = data{i, 'session'}{1};
    unit = data{i, 'unit'};
    
    try
        % fit models
        filename = fullfile(getenv('SSD'), 'paper2', 'modelling', 'glms', 'residual_glms', ...
            [session '_cell_' num2str(unit) '_glm.mat']);
        if overwrite || ~exist(filename, 'file')
            fprintf('(%3i/%i) %s, cell %3i: fitting GLMs\n', i, height(data), session, unit);
            fitResidualGlm(session, unit, 'verbose', false, 'parallel', false, 'save', true);
        end

        % plot
        plotFilename = fullfile(getenv('SSD'), 'paper2', 'modelling', 'glms', 'residual_glms', ...
            [session '_cell_' num2str(unit) '_glm.png']);
        if (overwrite || ~exist(plotFilename, 'file')) && exist(filename, 'file')
            plotResidualGlms(session, unit, 'save', true, 'visible', false)
        end
    catch exception
        fprintf('%s: PROBLEM! -> %s, unit %i\n', session, unit, exception.identifier)
    end
end
fprintf('\nfinished in %.1f minutes\n', toc/60)

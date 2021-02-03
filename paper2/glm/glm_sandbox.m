%% play around with GLMs :)


%% train all residual GLMs

overwrite = true;

[sessions, neurons] = getEphysSessions();
sessions = repelem(sessions, cellfun(@length, neurons));
neurons = cat(1, neurons{:});

skipInds = [7 27 61];


tic; fprintf('\nfitting residual GLMs for %i neurons...\n', length(sessions))

parfor i = 1:length(sessions)  % individual sessions are repeated for each neuron in session
    if ~ismember(i, skipInds)
        try
            % fit models
            filename = fullfile(getenv('SSD'), 'paper2', 'modelling', 'glms', 'residual_glms', ...
                [sessions{i} '_cell_' num2str(neurons(i)) '_glm.mat']);
            if overwrite || ~exist(filename, 'file')
                fprintf('(%3i/%i) %s, cell %3i: fitting GLMs\n', i, length(sessions), sessions{i}, neurons(i));
                fitResidualGlm(sessions{i}, neurons(i), 'verbose', false, 'parallel', false, 'save', true);
            end

            % plot
            plotFilename = fullfile(getenv('SSD'), 'paper2', 'modelling', 'glms', 'residual_glms', ...
                [sessions{i} '_cell_' num2str(neurons(i)) '_glm.png']);
            if (overwrite || ~exist(plotFilename, 'file')) && exist(filename, 'file')
                plotResidualGlms(sessions{i}, neurons(i), 'save', true, 'visible', false)
            end
        catch exception
            fprintf('%s: PROBLEM! -> %s\n', sessions{i}, exception.identifier)
        end
    end
end
fprintf('\nfinished in %.1f minutes\n', toc/60)

%% UPPER_LOWER GLMS

%% fit single upper_lower glm
session = '200622_000'; unit = 264;
[models, fitdata] = fitNeuronGlm(session, unit, 'parallel', true, 'save', true);
filename = fullfile(getenv('SSD'), 'paper2', 'modelling', 'glms', 'upper_lower_glms', ...
    [session '_cell_' num2str(unit) '_glm']);
plotGlmModels(session, unit, 'outputFileName', filename, 'visible', false)

%% train upper_lower GLMs for all sessions

overwrite = false;

[sessions, neurons] = getEphysSessions();
sessions = repelem(sessions, cellfun(@length, neurons));
neurons = cat(1, neurons{:});

skipInds = [7 27 61];

tic; fprintf('\nfitting models for %i neurons...\n', length(sessions))
parfor i = 1:length(sessions)
    if ~ismember(i, skipInds)
        try
            filename = fullfile(getenv('SSD'), 'paper2', 'modelling', 'glms', 'upper_lower_glms', ...
                [sessions{i} '_cell_' num2str(neurons(i)) '_glm.mat']);
            if overwrite || ~exist(filename, 'file')
                % train model
                fprintf('(%3i/%i) %s, cell %3i: fitting GLMs\n', i, length(sessions), sessions{i}, neurons(i));
                fitNeuronGlm(sessions{i}, neurons(i), 'verbose', false, ...
                    'method', 'refit', 'parallel', false, 'outputFileName', filename);

                % plot
                filename = fullfile(getenv('SSD'), 'paper2', 'modelling', 'glms', 'upper_lower_glms', ...
                    [sessions{i} '_cell_' num2str(neurons(i)) '_glm']);
                plotGlmModels(sessions{i}, neurons(i), 'outputFileName', filename, 'visible', false)
            end

        catch exception
            fprintf('%s: PROBLEM! -> %s\n', sessions{i}, exception.identifier)
        end
    end
end
fprintf('\nfinished in %.1f minutes\n', toc/60)


%% compute session durations

for i = 1:length(neurons)
    session = neurons{i}{1};
    load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'designMatrices', [session '_designMatrix.mat']), ...
        'dmat', 't', 'reward_all')
    duration = (t(end) - t(1)) / 60;
    fprintf('%s: %.1f minutes\n', session, duration)
end


%% compute single model and show predictions

session = '181020_001'; 83; neuron = 66;

% make new design matrix (use this to play around with different group settings)
filename = fullfile(getenv('SSD'), 'paper2', 'modelling', 'designMatrices', [session '_designMatrix.mat']);
makeDesignMatrix(session, 'timeDegrees', 3, 'outputFileName', filename);

% fit (and save) model
[models, fitdata] = fitNeuronGlm(session, neuron, 'method', 'refit', 'parallel', true, 'outputFileName', ...
    fullfile(getenv('SSD'), 'paper2', 'modelling', 'glms', 'upper_lower_glms', ...
    [session '_cell_' num2str(neuron) '_glm.mat']));

plotGlmModels(session, neuron);


%% compare refit, shuffle, and mask importance analyses

session = '181020_001';
neuron = 66;
parallel = true;

% make new design matrix (use this to play around with different group settings)
filename = fullfile(getenv('SSD'), 'paper2', 'modelling', 'designMatrices', [session '_designMatrix.mat']);
makeDesignMatrix(session, 'timeDegrees', 3, 'outputFileName', filename);

tic; [models_refit, fitdata] = fitNeuronGlm(session, neuron, 'method', 'refit', 'parallel', parallel, 'outputFileName', ...
    fullfile(getenv('SSD'), 'paper2', 'modelling', 'glms', 'upper_lower_glms', ...
    [session '_cell_' num2str(neuron) '_glm.mat'])); toc
% tic; models_shuffle = fitNeuronGlm(session, neuron, 'method', 'shuffle', 'parallel', parallel, 'outputFileName', ...
%     fullfile(getenv('SSD'), 'paper2', 'modelling', 'glms', [session '_cell_' num2str(neuron) '_glm.mat'])); toc
% tic; models_mask = fitNeuronGlm(session, neuron, 'method', 'mask', 'parallel', parallel, 'outputFileName', ...
%     fullfile(getenv('SSD'), 'paper2', 'modelling', 'glms', [session '_cell_' num2str(neuron) '_glm.mat'])); toc

% plot
close all; figure('color', 'white', 'position', [2.00 1043.00-313 1278.00 313.00]); hold on
plot(fitdata.t, fitdata.yRate)
plot(fitdata.t, fitdata.yhat)
set(gca, 'xlim', mean(fitdata.t)+[-10 10])

figure('color', 'white', 'position', [2.00 1043.00 1278.00 313.00]);
ngroups = height(models_refit)-1;
models = {models_refit, models_shuffle, models_mask};
names = {'refit', 'shuffle', 'mask'};

for i = 1:3
    subplot(1,3,i); title(names{i}); hold on
    plot(repmat(1:ngroups+1,2,1), [models{i}.dev_in models{i}.dev_out]', 'color', [.6 .6 .6])
    slow = scatter(1:ngroups+1, models{i}.dev_out, 40, [1 .4 .4], 'filled');
    shigh = scatter(1:ngroups+1, models{i}.dev_in, 40, [.4 .4 1], 'filled');
    set(gca, 'XTick', 1:ngroups+1, 'ylim', [0 .5], ...
        'XTickLabel', ['full'; fitdata.groups], 'XTickLabelRotation', 40)
end
legend([shigh slow], 'high', 'low')

plotGlmModels(session, neuron)
set(gcf, 'position',[2.00 2.00 1278.00 645.00])

%% compare parallelization within or around fitHeuronGlm

session = '181020_001';
neurons = [37 54 65 66 69 83];
neurons = repmat(neurons,1,2);


% inner loop in parallelization
tic
for i = 1:(length(neurons))
    fitNeuronGlm(session, neurons(i), 'method', 'refit', 'parallel', true);
end
fprintf('inner loop time: %.2f minutes\n', toc/60);

% outer loop in parallelization
tic
parfor i = 1:length(neurons)
    fitNeuronGlm(session, neurons(i), 'method', 'refit', 'parallel', false);
end
fprintf('outer loop time: %.2f minutes\n', toc/60);

%% find best smoothing value for gross kinematic vars

% load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'predictors', [session '_predictors.mat']), 'predictors');
smoothings = linspace(0,1,100);


t = predictors.t{1};
dt = t(2) - t(1);
gross = zscore(predictors{'velocity_vel', 'data'}{1});
paw = zscore(predictors{'paw4RH_x', 'data'}{1});

maxCorrs = nan(1,length(smoothings));
for i = 1:length(smoothings)
    smoothed = smooth(gross, max(smoothings(i)/dt, 1));
    maxCorrs(i) = max(xcorr(smoothed, paw, round(1/dt), 'coeff'));
end


close all; figure('position', [2.00 722.00 1278.00 634.00]); hold on;
smoothing = .2;
smoothed = smooth(gross, max(1, smoothing/dt));
plot(t, smoothed);
plot(t, paw);
set(gca, 'xlim', 60*10+[4 14])

figure('position', [280.00 89.00 560.00 420.00]); plot(smoothings, maxCorrs)
corr(smoothed, paw')


%% fit epoch glms

%% fit single epoch glm
session = '200622_000'; unit = 264;
[models, fitdata] = fitEpochGlm(session, unit, 'parallel', true, 'save', true);


%% fit all epoch glms

data = getUnitInfo();

tic
for i = 1:height(data)
    try
        fitEpochGlm(data.session{i}, data.unit(i), 'parallel', false, 'closeFig', true);
    catch exception
        fprintf('%s (%i): PROBLEM! -> %s\n', data.session{i}, data.unit(i), exception.identifier)
    end
end
fprintf('\nfinished in %.1f minutes\n', toc/60)

%% plot reward vs predicted reward vs delta heatmaps

% settings
x = linspace(-2, 4, 200);
nclusters = 3;
colors = lines(nclusters);
pcutoff = .001;

% get data tbl and add reponse columns
data = getUnitInfo();
m = height(data);
n = length(x);
tbl = table(nan(m,n), nan(m,n), nan(m,n), nan(m,1), ...
    'VariableNames', {'resp', 'resp_pred', 'resp_delta', 'p'});
data = cat(2, data, tbl);


% get all responses
for i = 1:height(data)
    disp(i/height(data))
    
    modelname = sprintf('E:\\lab_files\\paper2\\modelling\\glms\\epoch_glms\\%s_cell_%i_glm.mat', data.session{i}, data.unit(i));
    
    if exist(modelname, 'file')
        d = load(modelname);
        r = load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'runAnalyzed', [data.session{i} '_runAnalyzed.mat']), ...
            'rewardTimes', 'omissionTimes');
        reward_all = [r.rewardTimes; r.omissionTimes];

        t    = d.fitdata.t;
        y    = d.fitdata.yRate;                    % actual firing rate
        yhat = d.models{'postreward', 'yhatout'};  % predicted from model trained outside rewards
        
        % normalize
        mn = mean(y);
        stdev = std(y);
        y = (y - mn) / stdev;
        yhat = (yhat - mn) / stdev;
        
        % remove prediction outliers
        yhat(yhat<min(y)) = min(y);
        yhat(yhat>max(y)) = max(y);

        % compute responses for each trial
        unit_resp  = interp1(t, y, reward_all + x);     % response matrix for unit
        unit_pred  = interp1(t, yhat, reward_all + x);  % predicted response matrix for unit
        unit_delta = unit_resp - unit_pred;
        
        % add mean responses to table
        data.resp(i,:) = nanmean(unit_resp, 1);
        data.resp_pred(i,:) = nanmean(unit_pred, 1);
        data.resp_delta(i,:) = nanmean(unit_delta, 1);

        % response significance
        pre  = nanmean(unit_resp(:,x<=0), 2);
        post = nanmean(unit_resp(:,x>0), 2);
        [~, data.p(i)] = ttest(pre, post);
    end
end

%%

close all;
figure('color', 'white', 'menubar', 'none', 'position', [2.00 722.00 1278.00 688.00])
bins = data.p < pcutoff;
nrows = sum(bins);
clims = [-5 5];
cmap = customcolormap([0 .5 1], [1 .2 .2; 1 1 1; .2 .2 1]);



% actual response
subplot(1,3,1)
imagesc(x, 1:nrows, data.resp(bins,:), clims)
colormap(cmap); colorbar

% predicted
subplot(1,3,2)
imagesc(x, 1:nrows, data.resp_pred(bins,:), clims)
colormap(cmap); colorbar

% actual - predicted
subplot(1,3,3)
imagesc(x, 1:nrows, data.resp_delta(bins,:), clims)
colormap(cmap); colorbar










































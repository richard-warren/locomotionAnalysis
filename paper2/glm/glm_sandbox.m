%% play around with GLMs :)

%% train glms for several neurons

neurons = {{'180917_002', 64}, ...   % pawKinematics
           {'181001_002', 93}, ...
           {'181020_001', 83}, ...   % grossKinematics
           {'200726_000', 269}, ...
           {'200708_000', 66}, ...   % reward
           {'200702_000', 269}, ...
           {'200703_000', 257}, ...  % whisker contact
           {'181020_001', 66}, ...
           {'200626_000', 262}, ...  % face
           {'200622_000', 21}, ...
           {'181103_000', 140}, ...  % obstacle
           {'181019_002', 67}, ...
           {'200622_000', 258}, ...  % vision
           {'200703_000', 275}, ...
           {'200725_000', 29}, ...
           };

tic; fprintf('\nfitting models for %i neurons...\n', length(neurons))
parfor i = 1:length(neurons)
    try
        fprintf('%s, cell %i\n', neurons{i}{1}, neurons{i}{2});
        
        % train model
        fitNeuronGlm(neurons{i}{1}, neurons{i}{2}, 'verbose', false, ...
            'method', 'refit', 'parallel', true, 'outputFileName', ...
            fullfile(getenv('SSD'), 'paper2', 'modelling', 'glms', [neurons{i}{1} '_cell_' num2str(neurons{i}{2}) '_glm.mat']));
        
        % plot
        filename = fullfile(getenv('SSD'), 'paper2', 'modelling', 'glms', [neurons{i}{1} '_cell_' num2str(neurons{i}{2}) '_glm']);
        plotGlmModels(neurons{i}{1}, neurons{i}{2}, 'outputFileName', filename, 'visible', false)
    
    catch exception
        fprintf('%s: PROBLEM! -> %s\n', neurons{i}{1}, exception.identifier)
    end
end
fprintf('\nfinished in %.1f minutes\n', toc/60)

%% plot glm results

parfor i = 1:length(neurons)
    try
        fprintf('plotting %s, cell %i\n', neurons{i}{1}, neurons{i}{2});
        filename = fullfile(getenv('SSD'), 'paper2', 'modelling', 'glms', [neurons{i}{1} '_cell_' num2str(neurons{i}{2}) '_glm']);
        plotGlmModels(neurons{i}{1}, neurons{i}{2}, 'outputFileName', filename, 'visible', false)
    catch exception
        fprintf('%s: PROBLEM! -> %s\n', neurons{i}{1}, exception.identifier)
    end
end




%% compute single model and show predictions

session = '181020_001'; neuron = 66;
[models, fitdata] = fitNeuronGlm(session, neuron, 'method', 'refit', 'parallel', true);

close all; figure('color', 'white', 'position', [2.00 1043.00 1278.00 313.00]); hold on
plot(fitdata.t, fitdata.yRate)
plot(fitdata.t, fitdata.yhat)


%% compare refit, shuffle, and mask importance analyses

session = '181020_001';
neuron = 66;

[models_refit, groups] = fitNeuronGlm(session, neuron, 'method', 'refit');
models_shuffle = fitNeuronGlm('181020_001', 66, 'method', 'shuffle');
models_mask = fitNeuronGlm('181020_001', 66, 'method', 'mask');

% plot
close all; figure('color', 'white', 'position', [2.00 1043.00 1278.00 313.00]);
ngroups = height(models_refit)-1;
models = {models_refit, models_shuffle, models_mask};
names = {'refit', 'shuffle', 'mask'};

for i = 1:3
    subplot(1,3,i); title(names{i}); hold on
    plot(repmat(1:ngroups+1,2,1), [models{i}.dev_in models{i}.dev_out]', 'color', [.6 .6 .6])
    slow = scatter(1:ngroups+1, models{i}.dev_out, 40, [1 .4 .4], 'filled');
    shigh = scatter(1:ngroups+1, models{i}.dev_in, 40, [.4 .4 1], 'filled');
    set(gca, 'XTick', 1:ngroups+1, 'XTickLabel', ['full'; groups], 'XTickLabelRotation', 40)
end
legend([shigh slow], 'high', 'low')

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











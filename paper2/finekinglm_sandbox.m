%% fit finekin glm
session = '181020_001'; neuron = 66;
models = fitFinekinGlm(session, neuron, 'parallel', true, 'save', true);;
plotFinekinGlms(session, neuron, 'save', 'true');

%% train all finekin GLMs

overwrite = false;

[sessions, neurons] = getEphysSessions();
sessions = repelem(sessions, cellfun(@length, neurons));
neurons = cat(1, neurons{:});

skipInds = [7 27 61];


tic; fprintf('\nfitting residual GLMs for %i neurons...\n', length(sessions))

parfor i = 1:length(sessions)  % individual sessions are repeated for each neuron in session
    if ~ismember(i, skipInds)
        try
            % fit models
            filename = fullfile(getenv('SSD'), 'paper2', 'modelling', 'glms', 'finekin_glms', ...
                [sessions{i} '_cell_' num2str(neurons(i)) '_glm.mat']);
            if overwrite || ~exist(filename, 'file')
                fprintf('(%3i/%i) %s, unit %3i: fitting GLMs\n', i, length(sessions), sessions{i}, neurons(i));
                fitFinekinGlm(sessions{i}, neurons(i), 'verbose', false, 'parallel', false, 'save', true);
                plotFinekinGlms(sessions{i}, neurons(i), 'save', true, 'visible', false)
            end
        catch exception
            fprintf('%s: PROBLEM! -> %s\n', sessions{i}, exception.identifier)
        end
    end
end
fprintf('\nfinished in %.1f minutes\n', toc/60)

%% load models

modelNames = {'gross', 'fine', 'grossfine', 'full'};
folder = 'E:\lab_files\paper2\modelling\glms\finekin_glms\';

data = getUnitNuclei;
tbl = table(nan(n,1), nan(n,1), nan(n,1), nan(n,1), 'VariableNames', modelNames);
data = cat(2, data, tbl);
n = height(data);

files = dir(fullfile(folder, '*.mat'));
files = {files.name};

for i = 1:length(files)
    model = load(fullfile(folder, files{i}));
    session = model.fitdata.session;
    unit = model.fitdata.neuron;
    dataRow = find(strcmp(data.session, session) & data.unit==unit);
    data{dataRow, modelNames} = model.models{modelNames, 'dev'}';
end

%% gross vs gross+fine scatters

% settings
xlims = [0 .4];
ylims = [0 .6];
axes = {'gross', 'grossfine'};

close all
figure('color', 'white', 'position', [2.00 1200.00 1278.00 210.00], 'menubar', 'none')

locations = unique(data.nucleus);
n = length(locations);
colors = lines(n);


for i = 1:n
    subplot(1,n,i); hold on
    bins = strcmp(data.nucleus, locations{i});
    gross = data.(axes{1})(bins);
    grossfine = data.(axes{2})(bins);
    scatter(gross, grossfine, [], colors(i,:), 'filled', ...
        'MarkerFaceAlpha', .4)
    plot(xlims, xlims, 'color', [0 0 0 .2])
    set(gca, 'xlim', xlims, 'ylim', ylims)
    
    title(locations{i})
    if i==1
        xlabel(axes{1})
        ylabel(axes{2})
    end
end

%% bar plots for gross, fine

























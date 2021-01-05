%% fit finekin glm
session = '200622_000'; neuron = 269;
% models = fitFinekinGlm(session, neuron, 'parallel', true, 'save', true);
tic; plotFinekinGlms(session, neuron, 'save', 'true'); toc

%% train all finekin GLMs

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
            filename = fullfile(getenv('SSD'), 'paper2', 'modelling', 'glms', 'finekin_glms', ...
                [sessions{i} '_cell_' num2str(neurons(i)) '_glm.mat']);
            if overwrite || ~exist(filename, 'file')
                fprintf('(%3i/%i) %s, unit %3i: fitting GLMs\n', i, length(sessions), sessions{i}, neurons(i));
%                 fitFinekinGlm(sessions{i}, neurons(i), 'verbose', false, 'parallel', false, 'save', true);
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


modelNamesAll = [modelNames cellfun(@(x) [x '_r2'], modelNames, 'UniformOutput', false)];
data = getUnitNuclei;
n = height(data);
inits = repelem({nan(n,1)}, length(modelNamesAll));
tbl = table(inits{:}, 'VariableNames', modelNamesAll);
data = cat(2, data, tbl);

files = dir(fullfile(folder, '*.mat'));
files = {files.name};

for i = 1:length(files)
    disp(i/length(files))
    model = load(fullfile(folder, files{i}));
    
    session = model.fitdata.session;
    unit = model.fitdata.neuron;
    dataRow = find(strcmp(data.session, session) & data.unit==unit);
    data{dataRow, modelNames} = model.models{modelNames, 'dev'}';
    data{dataRow, 'r2'} = r2;
    
    % compute r2
    dt = diff(model.fitdata.t(1:2));
    y = model.fitdata.yRate';
    for j = 1:length(modelNames)
        % compute r2 for full model
        yhat = exp(model.models{modelNames{j}, 'model'}{1}.fit_preval) / dt;
        data{dataRow, [modelNames{j} '_r2']} = corr(yhat, y)^2;
    end
end

%% gross vs gross+fine scatters

% settings
xlims = [0 .4];
ylims = [0 .6];
axes = {'gross', 'grossfine'};

% xlims = [0 .8];
% ylims = [0 .8];
% axes = {'gross_r2', 'grossfine_r2'};

figure('color', 'white', 'position', [2.00 1200.00 1278.00 210.00], 'menubar', 'none')

locations = unique(data.nucleus);
n = length(locations);
% colors = lines(n);


for i = 1:n
    subplot(1,n,i); hold on
    bins = strcmp(data.nucleus, locations{i});
    gross = data.(axes{1})(bins);
    grossfine = data.(axes{2})(bins);
    scatter(gross, grossfine, [], [0 0 0], 'filled', ...
        'MarkerFaceAlpha', .4)
    plot(xlims, xlims, 'color', [0 0 0 .2])
    set(gca, 'xlim', xlims, 'ylim', ylims)
    
    title(locations{i})
    if i==1
        xlabel(axes{1})
        ylabel(axes{2})
    end
end

%% show units that have improved with fine kinematics

close all
nunits = 20;
r2sort = true;  % whether to sort by r squared rather than deviance explain

if r2sort
    [~, sortInds] = sort(data.grossfine_r2-data.gross_r2, 'descend', 'MissingPlacement', 'last');
else  
    [~, sortInds] = sort(data.grossfine-data.gross, 'descend', 'MissingPlacement', 'last');
end



for i = 1:nunits
    session = data{sortInds(i),'session'}{1};
    unit = data{sortInds(i),'unit'};
    nucleus = data{sortInds(i),'nucleus'}{1};
    plotFinekinGlms(session, unit, 'save', false);
    xlabel(sprintf('%s unit %i (%s)', session, unit, nucleus), 'Interpreter', 'none')
end

%% test r^2 across all neurons

% settings
otherR2s = struct('witton', .257, 'chabrol', .17, 'churchland', .2);

% models = {'gross_r2', 'grossfine_r2', 'full_r2'};
% ylims = [0 .15];

models = {'gross', 'grossfine', 'full'};
ylims = [0 .3];

xlims = [0 .8];
binEdges = 0:.02:1;

close all
figure('color', 'white', 'position', [130.00 840.00 632.00 374.00], 'menubar', 'none'); hold on
colors = lines(length(models));

% hists = cell(1, length(models));
for i = 1:length(models)
    d = data{:, models{i}};
    histogram(d, binEdges, 'Normalization', 'probability', 'FaceColor', colors(i,:), 'EdgeColor', 'none');
    plot([1 1]*nanmean(d), [0 1], 'color', colors(i,:), 'LineWidth', 2);
end

% add other papers
% for i = fieldnames(otherR2s)'
%      x = otherR2s.(i{1});
%      plot([x x], ylims, '--', 'color', [0 0 0 .4], 'linewidth', 2);
%      text(x, ylims(2), i{1}, 'Rotation', 90, 'BackgroundColor', 'white', ...
%          'HorizontalAlignment', 'right', 'FontSize', 10)
% end

for i = 1:length(models); lns(i) = plot([nan nan], 'color', colors(i,:), 'LineWidth', 2); end % create dummy lines
legend(lns, models, 'Interpreter', 'none', 'box', 'off')
set(gca, 'ylim', ylims, 'xlim', xlims)
xlabel('r^2')


%% todo: show increase in accuracy vs time over session




















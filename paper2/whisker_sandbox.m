%% explore whisker responses


%% collect data for all sessions
% (for each session get whisker neural responses for all units and 
%  trialData, which has variables associated with each trial)

sessions = getEphysSessions();
n = length(sessions);

data = table(cell(n,1), cell(n,1), cell(n,1), cell(n,1), cell(n,1), ...
    'RowNames', sessions, 'VariableNames', {'responses', 'unit_ids', 'trialData', 'whiskerAngle', 'deviance_explained'});

fprintf('preparing session: ')
for i = 1:length(sessions)
    fprintf('%i ', i)

    % load unit responses 
    load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'responses', [sessions{i} '_responses.mat']), 'responses');
    responses = responses{'whiskerContact', 'response'}{1};

    % load metadata for each trial
    load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'expData', [sessions{i} '_expData.mat']), 'expData');
    trialData = struct2table(expData.data(1).sessions(1).trials);
    leftFirst = nan(height(trialData), 1);  % add column for whether left paw is leading
    for j = 1:height(trialData); leftFirst(j) = trialData{j, 'paws'}(2).isLeading; end
    trialData = cat(2, trialData, table(leftFirst, 'VariableNames', {'leftFirst'}));
    
    % load whisker angles
    load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'predictors', [sessions{i} '_predictors.mat']), 'predictors');
    whiskerAngle = predictors('whiskerAngle', {'data', 't'});
    
    % load deviance explained
    load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'neuralData', [sessions{i} '_neuralData.mat']), 'unit_ids');
    deviance_explained = nan(1, length(unit_ids));
    for j = 1:length(unit_ids)
        filename = fullfile(getenv('SSD'), 'paper2', 'modelling', 'glms', 'residual_glms', ...
            [sessions{i} '_cell_' num2str(unit_ids(j)) '_glm.mat']);
        if exist(filename, 'file'); load(filename, 'models'); end
        deviance_explained(j) = models{'whiskers', 'dev'};
    end
    
    % store    
    data{sessions{i}, 'responses'} = {responses};
    data{sessions{i}, 'unit_ids'} = {unit_ids};
    data{sessions{i}, 'trialData'} = {trialData};
    data{sessions{i}, 'whiskerAngle'} = {whiskerAngle};
    data{sessions{i}, 'deviance_explained'} = {deviance_explained};
end
unitsPerSession = cellfun(@(x) size(x,3), data.responses);

% get x time axis
load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'responses', [sessions{1} '_responses.mat']), 'responses');
xlims = responses{'whiskerContact', 'xLims'};  % assumes xlims same for all sessions
x = linspace(xlims(1), xlims(end), size(data{1, 'responses'}{1},2));

% get deviance explained sort inds
deviance_explained = cat(2, data.deviance_explained{:});
[~, sortInds] = sort(deviance_explained, 'descend');

fprintf('all done!\n')


%% binned responses

% settings
ncols = 24;
binVar = 'isLightOn';
nbins = 3;
colors = copper(nbins);
xlims = [-.1 .25];
ncolsWisk = 8;  % columns for whisker angle plots



% speed notes: addding pause btwn sessions slows down A LOT // hold on
% doesn't affect speed // sorting subplots doesn't affect speed // patch
% may even be faster than plotting error lines

% inits
close all
figure('name', ['whiskers_' binVar], 'color', 'white', 'position', [1.00 1.00 2560.00 1363.00]); hold on
nrows = ceil(sum(unitsPerSession) / ncols);

tic
for i = 1:length(sessions)
    disp(i/length(sessions))

    responses = data{i, 'responses'}{1};
    trialData = data{i, 'trialData'}{1};
    [~,~,bins] = histcounts(trialData.(binVar), nbins);
    
    % plot neural PSTHs
    for j = 1:size(responses, 3)
        plotInd = sum(unitsPerSession(1:i-1)) + j;
        subplot(nrows, ncols, find(plotInd==sortInds)); hold on  % sorted
        
        for k = 1:nbins
            resp = responses(bins==k,:,j);
            respMean = nanmean(resp,1);
            respStd  = nanstd(resp,1);
            patch([x fliplr(x)], [(-respStd+respMean) fliplr(respStd+respMean)], colors(k,:), ...
                'FaceAlpha', .25, 'EdgeColor', 'none')  % shaded error bars
            plot(x, respMean, 'color', colors(k,:), 'lineWidth', 2)  % mean
            set(gca, 'XLim', xlims, 'xcolor', 'none')
        end
        plot([0 0], ylim, 'color', [.2 .2 .2])  % vertical line at zero
        
        % add text
        title(sprintf('%s %i', ...
            sessions{i}, data{i, 'unit_ids'}{1}(j)), ...
            'interpreter', 'none', 'FontSize', 8, 'FontWeight', 'normal')
        ylims = ylim;
        text(xlims(2), ylims(2), sprintf('%.3f', data{i, 'deviance_explained'}{1}(j)), ...
            'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 8)
    end
end
set(gca, 'xcolor', get(gca, 'ycolor'))  % time axis visible only for final subplot
toc

saveas(gcf, fullfile(getenv('SSD'), 'paper2', 'plots', 'whiskers', ...
    ['whiskerResponses_binnedby_' binVar '.png']))


% whisker angle PSTHs
nrows = ceil(length(sessions) / ncolsWisk);
figure('name', ['whiskerAngles_' binVar], 'color', 'white', 'position', [2.00 722.00 1278.00 634.00]); hold on

for i = 1:length(sessions)
    
    responses = data{i, 'responses'}{1};
    trialData = data{i, 'trialData'}{1};
    [~,~,bins] = histcounts(trialData.(binVar), nbins);
    
    % get whisker angle PSTHs
    whiskerAngle = data{i, 'whiskerAngle'}{1}.data{1};
    whiskerAngle = whiskerAngle - nanmean(whiskerAngle);
    t = data{i, 'whiskerAngle'}{1}.t{1};
    events = data{i, 'trialData'}{1}.wiskContactTimes;  % times of whisker contacts
    eventBins = events+xlims(1)>t(1) & events+xlims(2)<t(end);  % only include events with windows falling inside range of t
    wiskResponses = interp1(t, whiskerAngle, x+events);
    
    % plot whisker PSTH
    subplot(nrows, ncolsWisk, i); hold on
    for k = 1:nbins
        resp = wiskResponses(bins==k,:);
        respMean = nanmean(resp,1);
        respStd  = nanstd(resp,1);
        plot([0 0], ylim, 'color', [.2 .2 .2])  % vertical line at zero
        patch([x fliplr(x)], [(-respStd+respMean) fliplr(respStd+respMean)], colors(k,:), ...
            'FaceAlpha', .2, 'EdgeColor', 'none')
        plot(x, respMean, 'color', colors(k,:), 'lineWidth', 2)
        set(gca, 'XLim', xlims, 'xcolor', 'none')
    end
end
set(gca, 'xcolor', get(gca, 'ycolor'))  % time axis visible only for final subplot

saveas(gcf, fullfile(getenv('SSD'), 'paper2', 'plots', 'whiskers', ...
    ['whiskerAngleResponses_binnedby_' binVar '.png']))

%% todo: sort by something! (deviance explained, MI...)




















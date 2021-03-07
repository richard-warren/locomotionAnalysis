%% assess tuning steps over obstacle...



%% get step tuning (via residual model // use to write vids on lab comp sorted by step importance)

data = getUnitInfo(true);
data = data(ismember(data.nucleus, {'fastigial', 'interpositus', 'dentate'}), :);

% add importance for every unit (based on additional deviance explained in residual glm models)
data.stepImportance = nan(height(data), 1);
for i = 1:height(data)
    disp(i/height(data))
    fname = fullfile('E:\lab_files\paper2\modelling\glms\residual_glms', ...
        sprintf('%s_cell_%i_glm.mat', data.session{i}, data.unit(i)));
    if exist(fname, 'file')
        load(fname, 'models');
        data.stepImportance(i) = max(0, models{'pawKinematics', 'dev'});
    end
end

save('Y:\loco\obstacleData\data_transfer\to_remote\dataWithStepImportance.mat', 'data')


%% get neural responses
data = getStepResponses();
save(fullfile(getenv('SSD'), 'paper2', 'modelling', 'stepData', 'stepData.mat'), 'data');

%% plot that ish

% settings
ncols = 24;  % divisible by length(paws) ideally
nrows = 14;
paws = 4;
binvar = 'vel';
nbins = 3;
colors = copper(nbins);
prctileBins = true;  % whether to bin by percentiles
stepOverOnly = true;  % whether to restrict to step over obs
plotPredicted = false;


% overwrite some settings for logical vars
if islogical(data{1,'stepData'}{1}.(binvar))
    colors = [.2 .2 .2; lines(1)];
    nbins = 2;
end

% todo: only include highest obstacles... // leading vs lagging // pre sated trials only // sort by deviance explained
% check cells gradation in the step over vs normal steps condition...


% inits
close all
fname = sprintf('psths%i_%s.png', 1, binvar);
figure('name', fname, 'color', 'white', 'position', [1.00 1.00 2560.00 1363.00]); hold on
labels = {'LH', 'LF', 'RF', 'RH'};
offset = 0; fignum = 1; % subplot index offset (used when making multiple figures)
unitsPerSession = cellfun(@(x) size(x,1), data.responses(:,1));
sessions = data.Properties.RowNames;
nx = size(data{1, 'responses'}{1}, 3);  % number of points in x axis
x = linspace(0, 1, nx);
unitColors = lines(sum(unitsPerSession));


tic
for i = 1:length(sessions)  % for sessions
    disp(i/length(sessions))
    nunits = unitsPerSession(i);
%     stepData = data{i, 'stepData'};
    
    for j = 1:nunits  % for units
        unitStartPlot = (sum(unitsPerSession(1:i-1))+j-1) * length(paws);  % start plot for this session
        
        for k = 1:length(paws)  % for paw
            plotInd = unitStartPlot + k - offset;
                
            % make new figure if necessary
            if plotInd > ncols*nrows
                offset = offset + ncols*nrows;
                plotInd = plotInd - ncols*nrows;
                fignum = fignum + 1;
                saveas(gcf, ['E:\lab_files\paper2\plots\steps\' fname])
                fname = sprintf('psths%i_%s.png', fignum, binvar);
                figure('name', fname(1:end-4), 'color', 'white', 'position', [1.00 1.00 2560.00 1363.00]); hold on
            end

            subplot(nrows, ncols, plotInd); hold on

            responses = data{i, 'responses'}{paws(k)};  % (unit X step X 'time')
            predicted = data{i, 'predicted'}{paws(k)};
    
            if ~isempty(responses)
                var = data{i,'stepData'}{paws(k)}.(binvar);  % var for binning responses
                
                % restrict to step over trials
                if stepOverOnly
                    isStepOver = data{i,'stepData'}{paws(k)}.isStepOver;
                    var = var(isStepOver);
                    responses = responses(:, isStepOver, :);
                    predicted = predicted(:, isStepOver, :);
                end
                
                % determine trial bins
                if islogical(var)
                    bins = var+1;
                elseif prctileBins
                    binEdgesPrctile = linspace(0, 100, nbins+1);  % percentile bin edges
                    binEdges = prctile(var, binEdgesPrctile);
                    [~,~,bins] = histcounts(var, binEdges);
                else
                    [~,~,bins] = histcounts(var, nbins);
                end
                
                % plot, omg
                for m = 1:nbins
                    resp = squeeze(responses(j,bins==m,:));
                    pred = squeeze(predicted(j,bins==m,:));
                    if size(resp,2)==1; resp=resp'; pred=pred'; end  % correct orientation if only one trial in this bin
                    n = sum(~isnan(resp(:,1)));
                    
                    % plot response
                    mean = nanmean(resp, 1);
                    if n>1
                        stdev = nanstd(resp, 1);
                        patch([x fliplr(x)], [(-stdev+mean) fliplr(stdev+mean)], colors(m,:), ...
                            'FaceAlpha', .25, 'EdgeColor', 'none')  % shaded error bars
                    end
                    plot(x, mean, 'color', colors(m,:), 'lineWidth', 2)  % mean
                    
                    % plot predicted
                    if plotPredicted
                        mean = nanmean(pred, 1);
                        plot(x, mean, '--', 'color', colors(m,:), 'lineWidth', 1)  % mean
                    end
                    
                    plot([0 0], ylim, 'color', [.2 .2 .2])  % vertical line at zero
                end
            end
            
            unitind = sum(unitsPerSession(1:i-1))+j;
            title(sprintf('%s %i', ...
                sessions{i}, data{i, 'unit_ids'}{1}(j)), ...
                'interpreter', 'none', 'FontSize', 8, 'FontWeight', 'normal', 'Color', unitColors(unitind,:))
            yticks = get(gca, 'ytick');
            set(gca, 'XLim', [0 1], 'XTick', [], 'YTick', yticks([1,end]))
            xlabel(sprintf('%s', labels{paws(k)}))
        end
    end
end
set(gca, 'xcolor', get(gca, 'ycolor'))  % time axis visible only for final subplot
saveas(gcf, ['E:\lab_files\paper2\plots\steps\' fname])  % save final figure
toc




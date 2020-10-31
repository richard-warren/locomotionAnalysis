function plotResidualGlms(session, neuron, varargin)
% plots true and predicted firing rate for full glm, plus glm with only
% each predictor group included... also shows upper and lower bound
% deviance explained for each predictor group


% settings
s.save = true;                               % whether to save output to file
s.visible = true;                            % whether figure is visible
colors = [.4 .4 .4; .4 .4 1; 1 .4 .4];       % true firing rate, upper bound, lower bound
s.outputFileName = fullfile(...
    getenv('SSD'), 'paper2', 'modelling', 'glms', 'residual_glms', ...
    [session '_cell_' num2str(neuron) '_glm']);

% inits
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs
events = {'reward_all', 'whiskerContact', 'lick'};
eventColors = lines(length(events));
folder = fullfile(getenv('SSD'), 'paper2', 'modelling');
load(fullfile(folder, 'glms', 'residual_glms', [session '_cell_' num2str(neuron) '_glm.mat']), 'models', 'fitdata');
load(fullfile(folder, 'predictors', [session '_predictors.mat']), 'predictors');
groups = models.Properties.RowNames;
ngroups = length(groups);
t = fitdata.t;
y = fitdata.yRate;
dt = t(2) - t(1);
r = predictors{'reward_all', 'data'}{1};  % reward times
r = r(r>t(1) & r<t(end));
xlims = r(floor(length(r)/2 + [0 1])) + [-1; 1];  % spanning two rewards



fig = figure('name', sprintf('%s, neuron %i', session, neuron), 'color', 'white', 'position', [2.00 2.00 1278.00 1354.00]); hold on
axes = cell(1,ngroups);
ylims = [0 prctile(y, 99.9)];


% importance
subplot(ngroups+1, 1, 1); hold on
plot(repmat(1:ngroups,2,1), [zeros(height(models),1) models{:,'dev'}]', 'color', colors(1,:), 'LineWidth', 1)
scatter(1:ngroups, models{:,'dev'}, 40, colors(3,:), 'filled');
set(gca, 'XTick', 1:ngroups, 'XTickLabel', groups, 'XTickLabelRotation', 20);
ylabel('deviance explained')
pos = get(gca, 'position');
pos(1) = .35; pos(3) = .3;
set(gca, 'position', pos, 'xlim', [0 ngroups+1], 'ylim', [0 .3], 'TickDir', 'out');

for j = find(~cellfun(@isempty, models{:,'model'}))'  % only groups for which models were trained
    yhat_pre = exp(models{j, 'model_pre'}{1}.fit_preval) / dt;
    yhat = exp(models{j, 'model'}{1}.fit_preval) / dt;
    
    % fill in nans for models that weren't trained on full data
    if ~all(models{j, 'bins'})
        [yhatTemp, yhat_preTemp] = deal(nan(size(t)));
        yhatTemp(models{j, 'bins'}) = yhat;
        yhat_preTemp(models{j, 'bins'}) = yhat_pre;
        yhat = yhatTemp;
        yhat_pre = yhat_preTemp;
    end
    
    % firing rate predictions
    axes{j} = subplot(ngroups+1, 1, j+1); hold on
    plot(t, y, 'color', colors(1,:))
    plot(t, yhat_pre, 'color', colors(2,:))
    plot(t, yhat, 'color', colors(3,:), 'LineWidth', 1.5);
    
    % plot events
    for k = 1:length(events)
        rewards = predictors{events{k}, 'data'}{1};
        scatter(rewards, repelem(ylims(2), length(rewards)), 10, eventColors(k,:), 'filled')
        plot(repmat(rewards',2,1), repmat(ylims',1,length(rewards)), 'color', eventColors(k,:))
    end
  
    % prettify
    set(gca, 'box', 'off', 'xlim', xlims, 'TickDir', 'out', 'ylim', ylims)
    if j<ngroups; set(gca, 'xcolor', 'none'); end
    ylabel(groups{j})
end
linkaxes([axes{:}])

if s.save
%     savefig([s.outputFileName '.fig']);
    saveas(gcf, [s.outputFileName '.png']);
end

if ~s.visible; close(fig); end


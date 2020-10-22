function plotGlmModels(session, neuron, varargin)
% plots true and predicted firing rate for full glm, plus glm with only
% each predictor group included... also shows upper and lower bound
% deviance explained for each predictor group


% settings
s.outputFileName = '';                       % whether to save output to file
s.visible = true;                            % whether figure is visible
colors = [.4 .4 .4; .4 .4 1; 1 .4 .4];       % true firing rate, upper bound, lower bound


% inits
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs
events = {'reward_all', 'whiskerContact', 'lick'};
eventColors = lines(length(events));
folder = fullfile(getenv('SSD'), 'paper2', 'modelling');
load(fullfile(folder, 'glms', [session '_cell_' num2str(neuron) '_glm.mat']), 'models', 'fitdata');
load(fullfile(folder, 'predictors', [session '_predictors.mat']), 'predictors');
groups = models.Properties.RowNames;
ngroups = length(groups);
t = fitdata.t;
y = fitdata.yRate;
dt = t(2) - t(1);
r = predictors{'reward_all', 'data'}{1};  % reward times
xlims = r(floor(length(r)/2 + [0 1])) + [-2; 2];  % spanning two rewards
if ~s.visible; vis = 'off'; else; vis = 'on'; end


fig = figure('color', 'white', 'position', [2.00 2.00 1278.00 1354.00]); hold on
cols = 5;
axes = cell(1,ngroups);
ylims = [0 prctile(y, 99.9)];
for j = 1:ngroups
    yhat_upper = exp(models{groups{j}, 'model_in'}{1}.fit_preval) / dt;
%     if j>1; yhat_lower = exp(models{groups{j}, 'model_out'}{1}.fit_preval) / dt; end
    
    % importance
    subplot(ngroups,cols,(j-1)*cols+1); hold on
    plot(repmat(1:ngroups,2,1), [models{:,'dev_in'} models{:,'dev_out'}]', 'color', colors(1,:), 'LineWidth', 1)
    if j==1; inds = 1:ngroups; else; inds = j; end
    scat_upper = scatter(inds, models{:,'dev_in'}(inds), 40, colors(2,:), 'filled');
    scat_lower = scatter(inds, models{:,'dev_out'}(inds), 40, colors(3,:), 'filled');
    if j==1; set(gca, 'XTick', 1:ngroups, 'XTickLabel', groups, 'XTickLabelRotation', 40); end
    ylabel('deviance explained')
    
    % firing rate predictions
    axes{j} = subplot(ngroups,cols,(j-1)*cols+[2:cols]); hold on
    props = {};
    plot(t, y, 'color', colors(1,:), props{:})
    plot(t, yhat_upper, 'color', colors(2,:), props{:})
%     plot(t, yhat_lower, 'color', colors(3,:), props{:})
    
    % plot events
%     ylims = ylim;
    for k = 1:length(events)
        rewards = predictors{events{k}, 'data'}{1};
        scatter(rewards, repelem(ylims(2), length(rewards)), 10, eventColors(k,:), 'filled')
        plot(repmat(rewards',2,1), repmat(ylims',1,length(rewards)), 'color', eventColors(k,:))
    end
  
    % prettify
    set(gca, 'box', 'off', 'xlim', xlims, 'TickDir', 'out', 'ylim', ylims)
    if j<ngroups; set(gca, 'xcolor', 'none'); end
    if j==ngroups; legend([scat_upper, scat_lower], 'upper bound', 'lower bound'); end
    ylabel(groups{j})
end
linkaxes([axes{:}])

if ~isempty(s.outputFileName)
    savefig([s.outputFileName '.fig']);
    saveas(gcf, [s.outputFileName '.png']);
end

if ~s.visible; close(fig); end


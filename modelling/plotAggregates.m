function plotAggregates(varargin)

% runs plotAggregate for all predictors! mwahahahahaha!

% settings
s.folder = fullfile(getenv('OBSDATADIR'), 'figures', 'modelling', 'aggregates');
s.miMin = .05;  % discard units with mutual information < s.miMin
s.mahaMax = 8;  % units with post prob < posteriorCutoff are not assigned a group!
s.zscoreRows = false;  % where to each cell response independently
s.hideUnclustered = false;  % whether to hide units that are > s.mahaMax away from their supposed group
s.visible = false;  % whether figure is visible


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end % reassign settings passed in varargin
load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'modelling', 'aggregates.mat'), 'aggregates')
if s.visible; vis = 'on'; else; vis = 'off'; end


for i = 1:height(aggregates)
    try
        fig = plotAggregate(aggregates(i,:), ...
            'miMin', s.miMin, 'mahaMax', s.mahaMax, 'zscoreRows', s.zscoreRows, 'hideUnclustered', s.hideUnclustered);
        saveas(fig, fullfile(s.folder, [aggregates.Properties.RowNames{i} '_aggregate.png']))
    catch
        fprintf('WARNING! problem with %s\n', aggregates.Properties.RowNames{i})
    end
end
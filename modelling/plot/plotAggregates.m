function plotAggregates(varargin)

% runs plotAggregate for all predictors! mwahahahahaha!

% settings
s.folder = fullfile(getenv('OBSDATADIR'), 'figures', 'modelling', 'aggregates');  % figure will be saved in this folder
s.miMin = .05;  % discard units with mutual information < s.miMin
s.mahaMax = 8;  % units with post prob < posteriorCutoff are not assigned a group!
s.zscoreRows = false;  % where to each cell response independently
s.hideUnclustered = false;  % whether to hide units that are > s.mahaMax away from their supposed group
s.visible = false;  % whether figure is visible
s.noPawGroups = true;  % whether to NOT group paw variables


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end % reassign settings passed in varargin
load(fullfile(getenv('SSD'), 'modelling', 'aggregates', 'aggregates.mat'), 'aggregates');

for i = 1:height(aggregates)
    predictor = aggregates.Properties.RowNames{i};
    
    % set predictor-specific settings
    if contains(predictor, 'paw') && s.noPawGroups
        mahaMaxTemp = 1000;
        nGroups = 1;
    else
        mahaMaxTemp = s.mahaMax;
        nGroups = [];
    end
    
    try
        fprintf('plotting aggregate responses for: %s\n', predictor)
        fig = plotAggregate(aggregates(i,:), ...
            'miMin', s.miMin, 'mahaMax', mahaMaxTemp, 'zscoreRows', s.zscoreRows, 'nGroups', nGroups, ...
            'hideUnclustered', s.hideUnclustered, 'visible', s.visible, 'suppressWarning', true);
        saveas(fig, fullfile(s.folder, [aggregates.Properties.RowNames{i} '_aggregate.png']))
    catch
        fprintf('WARNING! problem with %s\n', predictor)
    end
end
disp('all done!')


function getFeatureImportance(session, varargin)

% computes the importance of each predictor for each unit (e.g. mutual
% information, r^2) // does this by loading predictors.mat and
% neuralResponses.mat // saves results in 'modelling\importance.mat' within
% the session folder


% settings
s.contSmps = 10000;  % randomly sample s.contSmps of continuous data for mutual information


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end % reassign settings passed in varargin
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'modelling', 'predictors.mat'), 'predictors');
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'modelling', 'responses.mat'), 'responses');
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'neuralData.mat'), ...
    'spkRates', 'timeStamps', 'unit_ids');

for i = 1:length(unit_ids)
    
    fprintf('%s, cell %i: computing feature importance...\n', session, unit_ids(i))
    
    % check that 'predictors' and 'responses' have same rows
    if ~isequal(predictors.Properties.RowNames, responses(i).responses.Properties.RowNames)
        disp('WARNING! predictors and responses have different row names!')
        break;
    end
    
    % initialize table
    nRows = height(predictors);
    cellImportance = table(nan(nRows,1), nan(nRows,1), predictors.type, ...
        'VariableNames', {'mi', 'r2', 'type'}, 'RowNames', predictors.Properties.RowNames);
    
    for j = 1:height(predictors)
        
        if predictors.type(j)=='event'
            response = responses(i).responses.response{j};
            xLims = responses(i).responses.xLims(j,:);
            x = linspace(xLims(1), xLims(2), size(response,2));
            x = repmat(x, 1, size(response,1));
            y = reshape(response', [], 1);
            cellImportance.mi(j) = computeMI(x, y);
            
        elseif predictors.type(j)=='epoch'
            response = responses(i).responses.response{j};
            xLims = responses(i).responses.xLims(j,:);
            x = linspace(xLims(1), xLims(2), size(response,2));
            x = repmat(x, 1, size(response,1));
            y = reshape(response', [], 1);
            cellImportance.mi(j) = computeMI(x, y);
            
        elseif predictors.type(j)=='continuous'
            spkRate = interp1(timeStamps, spkRates(i,:), predictors.t{j});
            cellImportance.mi(j) = computeMI(spkRate, predictors.data{j});
        end
    end
    
    % save cell importances
    importance(i).cell = unit_ids(i);
    importance(i).importance = cellImportance;
end

save(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'modelling', 'importance.mat'), 'importance')
disp('all done!')



% ---------
% FUNCTIONS
% ---------

function mi = computeMI(x, y)
    % compute mutual information between continuous variables x and y
    
    xy = rmmissing([x(:), y(:)]);
    if ~isempty(xy)
        xy = datasample(xy, min(size(xy,1), s.contSmps), 'Replace', false);
        x = py.numpy.array(xy(:,1)).reshape(int8(-1), int8(1));
        y = xy(:,2);
        mi = double(py.sklearn.feature_selection.mutual_info_regression(x,y));
    else
        mi = nan;
    end
end

end




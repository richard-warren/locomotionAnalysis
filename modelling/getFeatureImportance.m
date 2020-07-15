function getFeatureImportance(session, varargin)

% computes the importance of each predictor for each unit (e.g. mutual
% information, r^2) // does this by loading predictors.mat and
% neuralResponses.mat // saves results in 'modelling\importance.mat' within
% the session folder // note that every row contains info FOR EVERY CELL,
% with the last dimension corresponding to the cell dimension


% settings
s.contSmps = 5000;  % randomly sample s.contSmps of continuous data for mutual information
s.eventLims = [-.25 .5];  % (s) for mutual information only compute within these limits
s.epochLims = [-.25 1.25];  % (fraction of epoch) for mutual information only compute within these limits


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end % reassign settings passed in varargin
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'modelling', 'predictors.mat'), 'predictors');
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'modelling', 'responses.mat'), 'responses');
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'neuralData.mat'), ...
    'spkRates', 'timeStamps', 'unit_ids');

% make table
nRows = height(predictors);
uUnits = length(unit_ids);
importance = table(nan(nRows,uUnits), nan(nRows,uUnits), nan(nRows,uUnits), predictors.type, predictors.include, ...
    'VariableNames', {'mi', 'r', 'p',  'type', 'include'}, 'RowNames', predictors.Properties.RowNames);


for i = 1:length(unit_ids)
    fprintf('%s, (%i/%i) cell %i: computing feature importance...\n', session, i, length(unit_ids), unit_ids(i))
    
    % check that 'predictors' and 'responses' have same rows
    if ~isequal(predictors.Properties.RowNames, responses.Properties.RowNames)
        disp('WARNING! predictors and responses have different row names!')
        return;
    end
    
    % initialize table
    for j = find(predictors.include)'
        
        if predictors.type(j)=='event'
            response = responses.response{j}(:,:,i);
            xLims = responses.xLims(j,:);
            x = linspace(xLims(1), xLims(2), size(response,2));
            
            % restrict x axis
            bins = x>=s.eventLims(1) & x<=s.eventLims(2);
            x = x(bins);
            response = rmmissing(response(:, bins));
            
            % mutual information
            if size(response,1) > 1
                xFlat = repmat(x, 1, size(response,1));
                yFlat = reshape(response', [], 1);
                importance.mi(j,i) = computeMI(xFlat, yFlat, true);  % if more than one trial, treat predictor as 'discrete' (see sklearn docs)
            end

            
        elseif predictors.type(j)=='epoch'
            response = responses.response{j}(:,:,i);
            xLims = responses.xLims(j,:);
            x = linspace(xLims(1), xLims(2), size(response,2));
            
            % restrict x axis
            bins = x>=s.epochLims(1) & x<=s.epochLims(2);
            x = x(bins);
            response = rmmissing(response(:, bins));
            
            % mutual information
            if size(response,1) > 1
                xFlat = repmat(x, 1, size(response,1));
                yFlat = reshape(response', [], 1);
                importance.mi(j,i) = computeMI(xFlat, yFlat, true);
            end
            
            
        elseif predictors.type(j)=='continuous'
            % mutual information
            spkRate = interp1(timeStamps, spkRates(i,:), predictors.t{j});
            importance.mi(j,i) = computeMI(spkRate, predictors.data{j}, false);
            
            % correlation
            temp = rmmissing([spkRate', predictors.data{j}']);
            importance.r(j,i) = corr(temp(:,1), temp(:,2));
        end
    end
end

save(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'modelling', 'importance.mat'), 'importance')
disp('all done!')



% ---------
% FUNCTIONS
% ---------

function mi = computeMI(x, y, discrete)
    % compute mutual information between continuous variables x and y //
    % 'discrete' is whether predictor is continuous or discrete
    
    xy = rmmissing([x(:), y(:)]);
    if ~isempty(xy)
        xy = datasample(xy, min(size(xy,1), s.contSmps), 'Replace', false);
        x = py.numpy.array(xy(:,1)).reshape(int8(-1), int8(1));
        y = xy(:,2);
        mi = double(py.sklearn.feature_selection.mutual_info_regression(x, y, discrete));
    else
        mi = nan;
    end
end

end




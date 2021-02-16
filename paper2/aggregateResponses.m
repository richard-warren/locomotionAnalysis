function aggregateResponses(varargin)

% aggregates responses to all predictors for all cells! every aggregate
% response is a matrix, where the rows are cells and the colums are the x
% axis (determined by predictor type) // also stores metadata for each
% cell in 'cellInfo', eg session and unit_number, as well as mutual
% information for each cell for each predictor

% x limits for continuous variables are the 'median' for the x limits for
% each predictor across sessions

% todo: speed up by doing all cells at once...

% settings
s.binNum = 100;  % x axis resolution
s.eventLims = [-.25 .5];
s.epochLims = [-.25 1.25];

% todo: add std, density... // check that everyone has the same
% predictors!!! // utilize 'include' field


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end % reassign settings passed in varargin
cellInfo = getUnitInfo();
sessions = unique(cellInfo.session, 'stable');
load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'responses', [sessions{1} '_responses.mat']), 'responses');  % load sample session
nrows = height(responses);  % todo: check that different responses don't have different heights!
nunits = height(cellInfo);
aggregates = table(cell(nrows,1), nan(nrows,2), repmat({nan(length(sessions),2)}, nrows, 1), responses.type, cell(nrows,1), ...
    'VariableNames', {'aggregate', 'xLims', 'xLimsAll', 'type', 'include'}, 'RowNames', responses.Properties.RowNames);


% find x limits for each continuous variable for each session
% (these will be used to compute global x limits for each predictor also compute total number of units)
fprintf('getting xLims and nunits for %i sessions: ', length(sessions))
for i = 1:length(sessions)
    fprintf('%i ', i)
    load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'responses', [sessions{i} '_responses.mat']), 'responses');
    for j = 1:height(responses)  % loop over predictors
        aggregates.xLimsAll{j}(i,:) = responses.xLims(j,:);
    end
end

% set x limits for all predictor types
for i = 1:height(aggregates)
    if aggregates.type(i)=='event'
        aggregates.xLims(i,:) = s.eventLims;
    elseif aggregates.type(i)=='epoch'
        aggregates.xLims(i,:) = s.epochLims;
    elseif aggregates.type(i)=='continuous'
        aggregates.xLims(i,:) = nanmedian(aggregates.xLimsAll{i},1);
    end
end
fprintf('\n')


% initializations
xEvent = linspace(s.eventLims(1), s.eventLims(2), s.binNum);
xEpoch = linspace(s.epochLims(1), s.epochLims(2), s.binNum);
aggregates.aggregate = repmat({nan(nunits, s.binNum)}, height(aggregates), 1);
aggregates.include = repmat({nan(nunits, 1)}, height(aggregates), 1);


previousSession = '';
for i = 1:height(cellInfo)  % for all units
    session = cellInfo.session{i};
    unit = cellInfo.unit(i);
    fprintf('%s, unit %i: (%i/%i) adding responses...\n', session, unit, i, height(cellInfo))
    
    % load session data if new session reached
    if ~strcmp(previousSession, session)
        load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'neuralData', [session '_neuralData.mat']), ...
            'spkRates', 'unit_ids');
        load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'responses', [session '_responses.mat']), 'responses');
        previousSession = session;
    end
    
    unitind = find(unit==unit_ids);
    unitMean = nanmean(spkRates(unitind,:));
    unitStd = nanstd(spkRates(unitind,:));

    % for each predictor
    for k = 1:nrows

        aggregates.include{k}(i) = responses.include(k);

        if responses.include(k)
            if responses.type(k)=='event'
                y = nanmean(responses.response{k}(:,:,unitind),1);
                y = (y-unitMean) / unitStd;
                if any(~isnan(y))
                    x = linspace(responses.xLims(k,1), responses.xLims(k,2), length(y));  % original x axis
                    aggregates.aggregate{k}(i,:) = interp1(x, y, xEvent);
                end

            elseif responses.type(k)=='epoch'
                y = nanmean(responses.response{k}(:,:,unitind),1);
                y = (y-unitMean) / unitStd;
                if any(~isnan(y))
                    x = linspace(responses.xLims(k,1), responses.xLims(k,2), length(y));  % original x axis
                    aggregates.aggregate{k}(i,:) = interp1(x, y, xEpoch);
                end

            elseif responses.type(k)=='continuous'
                y = responses.response{k}(:,unitind);
                y = (y-unitMean) / unitStd;
                if any(~isnan(y))
                    x = linspace(responses.xLims(k,1), responses.xLims(k,2), length(y));  % original x axis
                    xi = linspace(aggregates.xLims(k,1), aggregates.xLims(k,2), s.binNum);
                    validInds = ~isnan(y);
                    temp = interp1(x(validInds), y(validInds), xi, 'linear');
                    temp = fillmissing(temp, 'linear', 'EndValues', 'nearest');  % in case response doesn't cover domain, extend on left and right...
                    aggregates.aggregate{k}(i,:) = temp;
                end
            end
        end
    end
end


file = fullfile(getenv('SSD'), 'paper2', 'modelling', 'response_aggregates.mat');
fprintf('saving results to: %s\n', file);
save(file, 'aggregates', 'cellInfo');



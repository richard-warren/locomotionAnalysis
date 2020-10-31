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
s.sessions = {};

% todo: add std, density... // check that everyone has the same
% predictors!!! // utilize 'include' field


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end % reassign settings passed in varargin
if isempty(s.sessions); s.sessions = getEphysSessions(); end
load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'responses', [s.sessions{1} '_responses.mat']), 'responses');  % load sample session
nRows = height(responses);
aggregates = table(cell(nRows,1), nan(nRows,2), repmat({nan(length(s.sessions),2)}, nRows, 1), responses.type, cell(nRows,1), ...
    'VariableNames', {'aggregate', 'xLims', 'xLimsAll', 'type', 'include'}, 'RowNames', responses.Properties.RowNames);


% find x limits for each continuous variable for each session
% (these will be used to compute global x limits for each predictor also compute total number of units)
fprintf('getting xLims and nUnits for %i sessions: ', length(s.sessions))
nUnits = 0;
for i = 1:length(s.sessions)
    fprintf('%i ', i)
    load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'responses', [s.sessions{i} '_responses.mat']), 'responses');
    responseSize = size(responses.response{find(responses.include,1,'first')});
    
    nUnits = nUnits + responseSize(end);  % different units are stored in the last dimension of responses
    
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
aggregates.aggregate = repmat({nan(nUnits, s.binNum)}, height(aggregates), 1);
aggregates.include = repmat({nan(nUnits, 1)}, height(aggregates), 1);
cellInfo = table(cell(nUnits,1), nan(nUnits,1), ...
    'VariableNames', {'session', 'unit'});



% for all sessions
rowInd = 0;
for i = 1:length(s.sessions)
    
    fprintf('%s: (%i/%i) adding responses...\n', s.sessions{i}, i, length(s.sessions))
    load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'neuralData', [s.sessions{i} '_neuralData.mat']), ...
        'spkRates', 'unit_ids');
    load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'responses', [s.sessions{i} '_responses.mat']), 'responses');
	
    % for all cells
    for j = 1:length(unit_ids)
        rowInd = rowInd + 1;
        cellMean = nanmean(spkRates(j,:));
        cellStd = nanstd(spkRates(j,:));

        % for each predictor
        for k = 1:height(responses)
            
            aggregates.include{k}(rowInd) = responses.include(k);
            
            if responses.include(k)
                if responses.type(k)=='event'
                    y = nanmean(responses.response{k}(:,:,j),1);
                    y = (y-cellMean) / cellStd;
                    if any(~isnan(y))
                        x = linspace(responses.xLims(k,1), responses.xLims(k,2), length(y));  % original x axis
                        aggregates.aggregate{k}(rowInd,:) = interp1(x, y, xEvent);
                    end

                elseif responses.type(k)=='epoch'
                    y = nanmean(responses.response{k}(:,:,j),1);
                    y = (y-cellMean) / cellStd;
                    if any(~isnan(y))
                        x = linspace(responses.xLims(k,1), responses.xLims(k,2), length(y));  % original x axis
                        aggregates.aggregate{k}(rowInd,:) = interp1(x, y, xEpoch);
                    end

                elseif responses.type(k)=='continuous'
                    y = responses.response{k}(:,j);
                    y = (y-cellMean) / cellStd;
                    if any(~isnan(y))
                        x = linspace(responses.xLims(k,1), responses.xLims(k,2), length(y));  % original x axis
                        xi = linspace(aggregates.xLims(k,1), aggregates.xLims(k,2), s.binNum);
                        validInds = ~isnan(y);
                        temp = interp1(x(validInds), y(validInds), xi, 'linear');
                        temp = fillmissing(temp, 'linear', 'EndValues', 'nearest');  % in case response doesn't cover domain, extend on left and right...
                        aggregates.aggregate{k}(rowInd,:) = temp;
                    end
                end
            end
        end
        
    cellInfo.session{rowInd} = s.sessions{i};
    cellInfo.unit(rowInd) = unit_ids(j);
    end
end

file = fullfile(getenv('SSD'), 'paper2', 'modelling', 'response_aggregates.mat');
fprintf('saving results to: %s\n', file);
save(file, 'aggregates', 'cellInfo');



function aggregateResponses(varargin)

% aggregates responses to all predictors for all cells! every aggregate
% response is a matrix, where the rows are cells and the colums are the x
% axis (determined by predictor type) // also stores metadata for each
% cell in 'cellInfo', eg session and unit_number, as well as mutual
% information for each cell for each predictor

% x limits for continuous variables are the 'median' for the x limits for
% each predictor across sessions

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
load(fullfile(getenv('OBSDATADIR'), 'sessions', s.sessions{1}, 'modelling', 'predictors.mat'), 'predictors');  % load sample session
nRows = height(predictors);
aggregates = table(cell(nRows,1), cell(nRows,1), nan(nRows,2), ...
    repmat({nan(length(s.sessions),2)}, nRows, 1), predictors.type, cell(nRows,1), ...
    'VariableNames', {'aggregate', 'mi', 'xLims', 'xLimsAll', 'type', 'include'}, ...
    'RowNames', predictors.Properties.RowNames);


% find x limits for each continuous variable for each session
% (these will be used to compute global x limits for each predictor also compute total number of units)
fprintf('getting xLims and nUnits for %i sessions: ', length(s.sessions))
nUnits = 0;
for i = 1:length(s.sessions)
    fprintf('%i ', i)
    load(fullfile(getenv('OBSDATADIR'), 'sessions', s.sessions{i}, 'modelling', 'responses.mat'), 'responses');
    
    nUnits = nUnits + length(responses);
    
    for j = 1:height(responses(1).responses)  % loop over predictors
        aggregates.xLimsAll{j}(i,:) = responses(1).responses.xLims(j,:);
    end
end

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
aggregates.aggregate = repmat({nan(nUnits, s.binNum)}, height(predictors), 1);
aggregates.mi = repmat({nan(nUnits, 1)}, height(predictors), 1);
aggregates.include = repmat({nan(nUnits, 1)}, height(predictors), 1);
cellInfo = table(cell(nUnits,1), nan(nUnits,1), ...
    'VariableNames', {'session', 'unit'});
rowInd = 0;


% for all sessions
for i = 1:length(s.sessions)
    
    fprintf('%s: (%i/%i) adding responses...\n', s.sessions{i}, i, length(s.sessions))
    load(fullfile(getenv('OBSDATADIR'), 'sessions', s.sessions{i}, 'neuralData.mat'), 'spkRates');
    load(fullfile(getenv('OBSDATADIR'), 'sessions', s.sessions{i}, 'modelling', 'responses.mat'), 'responses');
    load(fullfile(getenv('OBSDATADIR'), 'sessions', s.sessions{i}, 'modelling', 'importance.mat'), 'importance');

    % for all cells
    for j = 1:length(responses)
        rowInd = rowInd + 1;
        cellResponses = responses(j).responses;
        cellMean = nanmean(spkRates(j,:));
        cellStd = nanstd(spkRates(j,:));

        % for each predictor
        for k = 1:height(cellResponses)
            
            aggregates.include{k}(rowInd) = cellResponses.include(k);
            
            if cellResponses.include(k)
                if cellResponses.type(k)=='event'
                    y = nanmean(cellResponses.response{k},1);
                    y = (y-cellMean) / cellStd;
                    if any(~isnan(y))
                        x = linspace(cellResponses.xLims(k,1), cellResponses.xLims(k,2), length(y));  % original x axis
                        aggregates.aggregate{k}(rowInd,:) = interp1(x, y, xEvent);
                    end

                elseif cellResponses.type(k)=='epoch'
                    y = nanmean(cellResponses.response{k},1);
                    y = (y-cellMean) / cellStd;
                    if any(~isnan(y))
                        x = linspace(cellResponses.xLims(k,1), cellResponses.xLims(k,2), length(y));  % original x axis
                        aggregates.aggregate{k}(rowInd,:) = interp1(x, y, xEpoch);
                    end

                elseif cellResponses.type(k)=='continuous'
                    y = cellResponses.response{k};
                    y = (y-cellMean) / cellStd;
                    if any(~isnan(y))
                        x = linspace(cellResponses.xLims(k,1), cellResponses.xLims(k,2), length(y));  % original x axis
                        xi = linspace(aggregates.xLims(k,1), aggregates.xLims(k,2), s.binNum);
                        validInds = ~isnan(y);
                        temp = interp1(x(validInds), y(validInds), xi, 'linear');
                        temp = fillmissing(temp, 'linear', 'EndValues', 'nearest');  % in case response doesn't cover domain, extend on left and right...
                        aggregates.aggregate{k}(rowInd,:) = temp;
                    end
                end
                aggregates.mi{k}(rowInd) = importance(j).importance.mi(k);
            end
        end
        
    cellInfo.session{rowInd} = s.sessions{i};
    cellInfo.unit(rowInd) = responses(j).cell;
    end
end

file = fullfile(getenv('OBSDATADIR'), 'matlabData', 'modelling', 'aggregates.mat');
fprintf('saving results to: %s\n', file);
save(file, 'aggregates', 'cellInfo');



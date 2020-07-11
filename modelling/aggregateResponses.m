% function aggregateResponses(varargin)

% for each predictor, want matrices: response, std, density
% and vectors: mutual information


s.binNum = 100;
s.eventLims = [-.25 .5];
s.epochLims = [-.25 1.25];

% todo: determine total number of cells in advance? // figure out xLims in
% advance or manually define per predictor... // each heatmap has table
% with metadata, eg MI, correlation, std, x confidence...
% todo: check that everyone has the same predictors!!!

% initialize table
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end % reassign settings passed in varargin
sessions = getEphysSessions();
load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{1}, 'modelling', 'predictors.mat'), 'predictors');  % load sample session
nRows = height(predictors);
aggregates = table(cell(nRows,1), cell(nRows,1), nan(nRows,2), ...
    repmat({nan(length(sessions),2)}, nRows, 1), predictors.type, ...
    'VariableNames', {'aggregate', 'mi', 'xLims', 'xLimsAll', 'type'}, ...
    'RowNames', predictors.Properties.RowNames);

% find x limits for each continuous variable for each session
% these will be used to compute global x limits for each predictor
% also compute total number of units
sessions = getEphysSessions();
% sessions = sessions(1:6);  % temp
nUnits = 0;
for i = 1:length(sessions)
    disp(i)
    load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'modelling', 'responses.mat'), 'responses');
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

%%
xEvent = linspace(s.eventLims(1), s.eventLims(2), s.binNum);
xEpoch = linspace(s.epochLims(1), s.epochLims(2), s.binNum);
aggregates.aggregate = repmat({nan(nUnits, s.binNum)}, height(predictors), 1);
aggregates.mi = repmat({nan(nUnits, 1)}, height(predictors), 1);
cellInfo = table(cell(nUnits,1), nan(nUnits,1), ...
    'VariableNames', {'session', 'unit'});
rowInd = 0;

% for all sessions
for i = 1:length(sessions)
    
    fprintf('%s: adding responses from session...\n', sessions{i})
    load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'neuralData.mat'), 'spkRates');
    load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'modelling', 'responses.mat'), 'responses');
    load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'modelling', 'importance.mat'), 'importance');

    % for all cells
    for j = 1:length(responses)
        rowInd = rowInd + 1;
        cellResponses = responses(j).responses;
        cellMean = nanmean(spkRates(j,:));
        cellStd = nanstd(spkRates(j,:));

        % for each predictor
        for k = 1:height(cellResponses) %% only those included...

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
                    temp = interp1(x, y, xi, 'linear', 'extrap');
                    if all(isnan(temp)); keyboard; end
                    aggregates.aggregate{k}(rowInd,:) = temp;
                end
            end
            
            aggregates.mi{k}(rowInd) = importance(j).importance.mi(k);
        end
        
    cellInfo.session{rowInd} = sessions{i};
    cellInfo.unit(rowInd) = responses(j).cell;
    end
end




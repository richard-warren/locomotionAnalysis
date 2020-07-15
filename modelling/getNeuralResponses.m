function getNeuralResponses(session, varargin)

% for each predictor in predictors.mat, save (trial X xaxis) matrix of
% responses to that predictor // responses are handled differently for
% events, epochs, and continuous vars // for events, we store simple PSTH
% firing rates // for epochs, each epoch (which varies in length) is
% stretched over a common x axis // for continuous vars, TBD!!!



% settings
s.eventLims = [-1 1];  % (s)
s.epochLims = [-.5 1.5];  % (s)
s.epochGridNum = 200;  % number of points in epoch x axis
s.percentileLims = [.1 99.9];  % limits for continuous variables
s.contGridNum = 200;    % number of points in continuous x axis
s.contWindowSz = .05;  % width of moving average window, expressed as fraction of x axis
s.maxEpochs = 1000;  % if more than s.maxEpochs epochs, only compute central s.maxEpochs


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end % reassign settings passed in varargin
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'neuralData.mat'), ...
    'unit_ids', 'spkRates', 'timeStamps');
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'modelling', 'predictors.mat'), 'predictors');
dt = median(diff(timeStamps));
xEvent = s.eventLims(1) : dt : s.eventLims(2);
xEpoch = linspace(s.epochLims(1), s.epochLims(2), s.epochGridNum);


% for each cell
for i = 1:length(unit_ids)
    fprintf('%s: unit (%i/%i) %i: computing neural responses...\n', session, i, length(unit_ids), unit_ids(i))
    
    nRows = height(predictors);
    cellResponses = table(cell(nRows,1), cell(nRows,1), cell(nRows,1), predictors.type, nan(nRows,2), predictors.include, ...
        'VariableNames', {'response', 'std', 'density', 'type', 'xLims', 'include'}, ...
        'RowNames', predictors.Properties.RowNames);
    
    tLims = timeStamps([find(~isnan(spkRates(i,:)),1,'first') find(~isnan(spkRates(i,:)),1,'last')]);
    spkBins = ~isnan(spkRates(i,:));
    
    
    for j = find(predictors.include)'
        
        if predictors.type(j)=='event'
            events = predictors.data{j};
            
            response = nan(length(events), length(xEvent));
            for k = 1:length(events)
                startInd = find(timeStamps>=(events(k)+s.eventLims(1)),1,'first');
                endInd = startInd+length(xEvent)-1;
                if ~isempty(startInd) && endInd<=size(spkRates,2)
                    response(k,:) = spkRates(i,startInd:startInd+length(xEvent)-1);
                end
            end
            response = fillmissing(response, 'linear', 'EndValues', 'nearest');  % i don't think it's possible for NaNs to appear here, but just to be safe...
            cellResponses.response{j} = response;
            cellResponses.xLims(j,:) = s.eventLims;
            if all(isnan(response(:))); cellResponses.include(j) = 0; end  % this will happen when no events occur during the unit's usable period

            
        elseif predictors.type(j)=='epoch'
            epochs = predictors.data{j};
            epochs = epochs(all(epochs>tLims(1) & epochs<tLims(2),2),:);
            response = nan(size(epochs,1), length(xEpoch));
            durations = diff(epochs,1,2);
            
            % if too many epochs find the elements of central duration
            if size(epochs,1)>s.maxEpochs
                middleIndsSorted = floor(length(durations)/2 - s.maxEpochs*.5) + (1:s.maxEpochs);
                [~, sortInds] = sort(durations);
                inds = sortInds(middleIndsSorted)';  % !!! should maybe pick random intead of central-duration epochs
            else
                inds = 1:length(epochs);
            end
            
            for k = inds
                epoch = durations(k)*s.epochLims + epochs(k,1);
                epochBins = timeStamps>=epoch(1) & timeStamps<=epoch(2) * spkBins;
                if any(epochBins)
                    response(k,:) = interp1(timeStamps(epochBins), spkRates(i,epochBins), ...
                        linspace(epoch(1), epoch(2), s.epochGridNum), 'linear');
                    response(k,:) = fillmissing(response(k,:), 'linear', 'EndValues', 'nearest');
                end
            end
            cellResponses.response{j} = response;
            cellResponses.xLims(j,:) = s.epochLims;
            if all(isnan(response(:))); cellResponses.include(j) = 0; end  % this will happen when no events occur during the unit's usable period
            
        
        elseif predictors.type(j)=='continuous'
            spkRate = interp1(timeStamps, spkRates(i,:), predictors.t{j});
            spkRateBins = ~isnan(spkRate);
            xLims = prctile(predictors.data{j}, s.percentileLims);
            xGrid = linspace(xLims(1), xLims(2), s.contGridNum);
            winSz = s.contWindowSz*diff(xLims);
            
            [response, responseStd, density] = deal(nan(1, length(xGrid)));
            for k = 1:length(xGrid)
                bins = predictors.data{j} >= (xGrid(k)-.5*winSz) & ...
                    predictors.data{j} < (xGrid(k)+.5*winSz) & ...
                    spkRateBins;
                if any(bins)
                    response(k) = sum(spkRate(bins)) / sum(bins);
                    responseStd(k) = std(spkRate(bins));
                    density(k) = sum(bins);  % !!! density should not be affected by spkRatebins
                end
            end
            response = fillmissing(response, 'linear', 'EndValues', 'nearest');
            if all(isnan(response)); cellResponses.include(j) = 0; end  % not sure if this should every happen really...
            density = density / sum(density);
            
            cellResponses.response{j} = response;
            cellResponses.xLims(j,:) = xLims;
            cellResponses.std{j} = responseStd;
            cellResponses.density{j} = density;
        end 
    end
    
    % save cell responses
    responses(i).cell = unit_ids(i);
    responses(i).responses = cellResponses;
end

save(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'modelling', 'responses.mat'), 'responses', '-v7.3')
disp('all done!')


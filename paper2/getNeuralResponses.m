function getNeuralResponses(session, varargin)

% for each predictor in predictors.mat, save (trial X x X cell) matrix of
% responses to that predictor // responses are handled differently for
% events, epochs, and continuous vars // for events, we store simple PSTH
% firing rates // for epochs, each epoch (which varies in length) is
% stretched over a common x axis // for continuous moving average // note 
% that every row contains info FOR EVERY CELL, with the last dimension 
% corresponding to the cell dimension



% settings
s.eventLims = [-1 1];  % (s)
s.epochLims = [-.5 1.5];  % (s)
s.gridNum = 200;  % number of points in epoch and event x axis
s.percentileLims = [1 99];  % limits for continuous variables
s.contGridNum = 200;    % number of points in continuous x axis
s.contWindowSz = .05;  % width of moving average window, expressed as fraction of x axis
s.maxEpochs = 500;  % if more than s.maxEpochs epochs, only compute central s.maxEpochs


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end % reassign settings passed in varargin
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'neuralData.mat'), ...
    'unit_ids', 'spkRates', 'timeStamps');
load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'predictors', [session '_predictors.mat']), 'predictors');

% initialize table
nRows = height(predictors);
responses = table(cell(nRows,1), cell(nRows,1), cell(nRows,1), predictors.type, nan(nRows,2), predictors.include, ...
    'VariableNames', {'response', 'std', 'density', 'type', 'xLims', 'include'}, ...
    'RowNames', predictors.Properties.RowNames);
t = predictors{find(predictors.type=='continuous',1,'first'), 't'}{1};  % (assumes same time grid for all predictors!)

% put spikes on same time grid as predictors
nUnits = length(unit_ids);
if nUnits==1
    spkRates = interp1(timeStamps, spkRates, t, 'linear');
else
    spkRates = interp2(timeStamps, [1:nUnits]', spkRates, t, [1:nUnits]', 'linear');
end
clear timeStamps

% define x axes
xEvent = linspace(s.eventLims(1), s.eventLims(2), s.gridNum);
xEpoch = linspace(s.epochLims(1), s.epochLims(2), s.gridNum);



fprintf('%s: getting neural responses...\n', session)
for i = find(predictors.include)'

    if predictors.type(i)=='event'
        events = predictors.data{i};
        response = nan(length(events), length(xEvent), nUnits);
        
        % temp
        eventBins = events+s.eventLims(1)>t(1) & events+s.eventLims(2)<t(end);  % only include events with windows falling inside range of t
        for j = 1:nUnits
            response(eventBins,:,j) = interp1(t, spkRates(j,:), xEvent+events(eventBins));  % compute firing rate for all trials in one shot
        end
        
        response = fillmissing(response, 'linear', 2, 'EndValues', 'nearest');  % i don't think it's possible for NaNs to appear here, but just to be safe...
        responses.response{i} = response;
        responses.xLims(i,:) = s.eventLims;
        % todo: determine if cell should be included here? as opposed to predictor?


    elseif predictors.type(i)=='epoch'
        epochs = predictors.data{i};
        response = nan(size(epochs,1), length(xEpoch), nUnits);
        durations = diff(epochs,1,2);

        % if too many epochs find the elements of central duration
        % (note) this may include epochs that are not within the usable
        % period for all units!
        if size(epochs,1)>s.maxEpochs
            middleIndsSorted = floor(length(durations)/2 - s.maxEpochs*.5) + (1:s.maxEpochs);
            [~, sortInds] = sort(durations);
            inds = sortInds(middleIndsSorted)';  % !!! should maybe pick random intead of central-duration epochs
        else
            inds = 1:length(epochs);
        end

        for j = inds
            epoch = durations(j)*s.epochLims + epochs(j,1);  % start and end times for this epoch
            epochBins = t>=epoch(1) & t<=epoch(2);
            if any(epochBins)
                if nUnits==1
                    response(j,:,:) = interp1(t(epochBins), spkRates(:,epochBins), ...
                        linspace(epoch(1), epoch(2), s.gridNum));
                else
                    response(j,:,:) = interp2(t(epochBins), ...
                        [1:nUnits]', spkRates(:,epochBins), ...
                        linspace(epoch(1), epoch(2), s.gridNum), [1:nUnits]', 'linear')';
                end
            end
        end
        response = fillmissing(response, 'linear', 2, 'EndValues', 'nearest');
        responses.response{i} = response;
        responses.xLims(i,:) = s.epochLims;
        % todo: determine if cell should be included here? as opposed to predictor?
        

    elseif predictors.type(i)=='continuous'
        xLims = prctile(predictors.data{i}, s.percentileLims);
        xGrid = linspace(xLims(1), xLims(2), s.contGridNum);
        winSz = s.contWindowSz*diff(xLims);

        [response, responseStd] = deal(nan(length(xGrid), nUnits));
        density = nan(1, length(xGrid));
        for j = 1:length(xGrid)  % could this be sped up with hist?
            xBins = predictors.data{i} >= (xGrid(j)-.5*winSz) & ...
                predictors.data{i} < (xGrid(j)+.5*winSz);
            density(j) = sum(xBins);
            if any(xBins)
                response(j,:) = nanmean(spkRates(:,xBins), 2);
                responseStd(j,:) = nanstd(spkRates(:,xBins), 0, 2);
            end
        end
        response = fillmissing(response, 'linear', 2, 'EndValues', 'nearest');  % shouldn't these never be necessary?
        responseStd = fillmissing(response, 'linear', 2, 'EndValues', 'nearest');
        density = density / sum(density);

        responses.response{i} = response;
        responses.xLims(i,:) = xLims;
        responses.std{i} = responseStd;
        responses.density{i} = density;
        % todo: determine if cell should be included here? as opposed to predictor?
    end 
end

save(fullfile(getenv('SSD'), 'paper2', 'modelling', 'responses', [session '_responses.mat']), 'responses', '-v7.3');
fprintf('%s: all done getting neural responses :)\n', session)


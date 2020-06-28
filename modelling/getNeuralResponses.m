function getNeuralResponses(session)

% for each predictor in predictors.mat, save (trial X xaxis) matrix of
% responses to that predictor // responses are handled differently for
% events, epochs, and continuous vars // for events, we store simple PSTH
% firing rates // for epochs, each epoch (which varies in length) is
% stretched over a common x axis // for continuous vars, TBD!!!

% todo: add duration to responses // add phase // add 'include' or
% something for features not successfully tracked?


% settings
s.eventLims = [-1 1];  % (s)
s.epochLims = [-.5 1.5];  % (s)
s.epochGridNum = 200;  % number of points in epoch x axis
s.percentileLims = [1 99];  % limits for continuous variables
s.contGridNum = 200;    % number of points in continuous x axis
s.contWindowSz = .05;  % width of moving average window, expressed as fraction of x axis
s.plotResponses = true;  % whether to make and save plot of all responses


% initializations
cellData = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'cellData.csv'));
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'neuralData.mat'), ...
    'unit_ids', 'spkRates', 'timeStamps');
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'modelling', 'predictors.mat'), 'predictors');
dt = median(diff(timeStamps));
xEvent = s.eventLims(1) : dt : s.eventLims(2);
xEpoch = linspace(s.epochLims(1), s.epochLims(2), s.epochGridNum);


% for each cell
for i = 1:length(unit_ids)
    fprintf('%s, cell %i: computing neural responses...\n', session, cellData.unit_id(i))
    
    cellResponses = table({}, categorical({}, {'event', 'epoch', 'continuous'}), [], ...
        'VariableNames', {'response', 'type', 'xLims'}, 'RowNames', {});
    tLims = timeStamps([find(~isnan(spkRates(1,:)),1,'first') find(~isnan(spkRates(1,:)),1,'last')]);
    
    if s.plotResponses
        figure('color', 'white', 'name', sprintf('%s_cell%i', session, cellData.unit_id(i)), 'position', [185.00 58.00 1583.00 915.00]);
        rows = ceil(sqrt(height(predictors)));  % same num row and cols
        yMax = prctile(spkRates(i,:), 99);
    end

    
    for j = 1:height(predictors)
        
        if predictors.type(j)=='event'
            events = predictors.data{j};
            try
            events = events(all((events+s.eventLims)>tLims(1) & (events+s.eventLims)<tLims(2),2));
            catch; keyboard; end
            
            response = nan(length(events), length(xEvent));
            for k = 1:length(events)
                startInd = find(timeStamps>=events(k)+s.eventLims(1),1,'first');
                response(k,:) = spkRates(i,startInd:startInd+length(xEvent)-1);
            end
            cellResponses = addResponse(cellResponses, ...
                predictors.Properties.RowNames{j}, response, 'event', s.eventLims);
            
            if s.plotResponses
                subplot(rows, rows, j); hold on
                plot([0 0], [0 yMax], 'color', [0 0 0 .4])
                respMean = nanmean(response,1);
                respStd = nanstd(response,1);
                plot(xEvent, respMean, 'LineWidth', 3, 'color', [.2 .2 1])
                plot(xEvent, respMean + [respStd; -respStd], 'LineWidth', 1, 'color', [.2 .2 1 .4])
                xlabel(predictors.Properties.RowNames{j}, 'Interpreter', 'none')
                set(gca, 'xlim', s.eventLims, 'ylim', [0 yMax])
                pause(.01)
            end
            
            
        elseif predictors.type(j)=='epoch'
            epochs = predictors.data{j};
            epochs = epochs(all(epochs>tLims(1) & epochs<tLims(2),2),:);
            response = nan(length(epochs), length(xEpoch));
            
            for k = 1:size(epochs,1)
                dt = (epochs(k,2)-epochs(k,1));
                epoch = dt*s.epochLims + epochs(k,1);
                epochBins = timeStamps>=epoch(1) & timeStamps<=epoch(2);
                response(k,:) = interp1(timeStamps(epochBins), spkRates(i,epochBins), ...
                    linspace(epoch(1), epoch(2), s.epochGridNum), 'linear', 'extrap');
            end
            cellResponses = addResponse(cellResponses, ...
                predictors.Properties.RowNames{j}, response, 'epoch', s.epochLims);
            
            if s.plotResponses
                subplot(rows, rows, j); hold on
                plot([0 0; 1 1]', [0 yMax; 0 yMax]', 'color', [0 0 0 .4])
                respMean = nanmean(response,1);
                respStd = nanstd(response,1);
                plot(xEpoch, respMean, 'LineWidth', 3, 'color', [1 .2 .2])
                plot(xEpoch, respMean + [respStd; -respStd], 'LineWidth', 1, 'color', [1 .2 .2 .4])
                xlabel(predictors.Properties.RowNames{j}, 'Interpreter', 'none')
                set(gca, 'xlim', s.epochLims, 'ylim', [0 yMax])
                pause(.01)
            end
                
        
        elseif predictors.type(j)=='continuous'
            spkRate = interp1(timeStamps, spkRates(i,:), predictors.t{j});
            xLims = prctile(predictors.data{j}, s.percentileLims);
            xGrid = linspace(xLims(1), xLims(2), s.contGridNum);
            winSz = s.contWindowSz*diff(xLims);
            
            response = nan(1, length(xGrid));
            for k = 1:length(xGrid)
                bins = predictors.data{j} >= (xGrid(k)-.5*winSz) & ...
                    predictors.data{j} < (xGrid(k)+.5*winSz) & ...
                    ~isnan(spkRate);
                if any(bins); response(k) = sum(spkRate(bins)) / sum(bins); end
            end
            cellResponses = addResponse(cellResponses, ...
                predictors.Properties.RowNames{j}, response, 'continuous', xLims);
            
            if s.plotResponses
                subplot(rows, rows, j); hold on
                inds = randperm(length(spkRate), 1000);
                scatter(predictors.data{j}(inds), spkRate(inds), 5, 'filled', 'markerfacealpha', .5)
                plot(xGrid, response, 'LineWidth', 3)
                xlabel(predictors.Properties.RowNames{j}, 'Interpreter', 'none')
                if xLims(1)~=xLims(2)  % somewhat of a hack to avoid constant predictors, which occur when feature is not successfully tracked...
                    set(gca, 'xlim', xLims, 'ylim', [0 yMax])
                end
                pause(.01)
            end
        end 
    end
    
    % save cell responses somehow...
    responses(i).cell = cellData.unit_id(i);
    responses(i).responses = cellResponses;
    if s.plotResponses
        saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures', 'modelling', 'responses', ...
            sprintf('%s cell%i responses.png', session, cellData.unit_id(i))));
    end
end

save(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'modelling', 'responses.mat'), 'responses')
disp('all done!')



function responses = addResponse(responses, name, data, type, xLims)
    % extend predictors table by adding new row
    
    if ~exist('xLims', 'var'); xLims = []; end
    newRow = table({data}, categorical({type}, {'event', 'epoch', 'continuous'}), xLims, ...
        'VariableNames', {'response', 'type', 'xLims'}, 'RowNames', {name});
    responses = [responses; newRow];
%     fprintf('adding: %s...\n', name)
end

end




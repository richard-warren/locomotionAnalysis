%% assess tuning to paw contacts...


%% collect data

% for session 
%   get 8 contact time types // model predictions // firing rates for all cells
%   for cell
%     psth for actual/predicted, for each of 8 contact types

% settings
minTouchFrames = 3;  % only include touches in contact with obstacle for this number of consecutive frames (this is actually an approximation only... i use a simple medfilt for this, which is close to but not the same as duration thresholding)
x = linspace(-.5,1,1000);  % x axis for PSTHs

sessions = getEphysSessions();
n = length(sessions);

cellInits = repmat({{}, {}, {}, {}}, n, 1);
data = table(cellInits, cellInits, cell(n,1), cell(n,1), 'RowNames', sessions, 'VariableNames', ...
    {'dorsal_touches', 'ventral_touches', 'responses', 'predicted'});

fprintf('preparing session: ')
for i = 1:length(sessions)
    fprintf('%i ', i)

    % load session touch times
    load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'runAnalyzed.mat'), ...
        'frameTimeStamps', 'touchesPerPaw', 'touchClassNames');
    dorsalInds = find(contains(touchClassNames, 'orsal'));
    ventralInds = find(contains(touchClassNames, 'entral'));
    dorsalTouches = ismember(touchesPerPaw, dorsalInds);
    ventralTouches = ismember(touchesPerPaw, ventralInds);
    
    % debounce
    dorsalTouches = logical(medfilt1(double(dorsalTouches), minTouchFrames*2-1));
    ventralTouches = logical(medfilt1(double(ventralTouches), minTouchFrames*2-1));
    
    % get doral and ventral touch times for each of four paws
    for j = 1:4
        data{i, 'dorsal_touches'}{j} = frameTimeStamps(find(diff(dorsalTouches(:,j))==1)+1);
        data{i, 'ventral_touches'}{j} = frameTimeStamps(find(diff(ventralTouches(:,j))==1)+1);
    end
    
    
    % get per cell PSTHs
    load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'neuralData', [sessions{i} '_neuralData.mat']), ...
        'timeStamps', 'spkRates', 'unit_ids');
    responses = cell(4, 2, length(unit_ids));  % (paw) X (dorsal ventral) X (neuron)
    predicted = cell(4, 2, length(unit_ids));  % (paw) X (dorsal ventral) X (neuron)
    
    for j = 1:length(unit_ids)
        % load model predictions
        file = fullfile(getenv('SSD'), 'paper2', 'modelling', 'glms', 'residual_glms', ...
            [sessions{i} '_cell_' num2str(unit_ids(j)) '_glm.mat']);
        hasYhat = exist(file, 'file');
        if hasYhat; load(file, 'fitdata'); else; fprintf('(no yhat) '); end
        
        % compute PSTHs (for actual and predicted firing rate)
        for k = 1:4
            
            dorsalEvents = data{i, 'dorsal_touches'}{k};
            if ~isempty(dorsalEvents)
                responses{k,1,j} = interp1(timeStamps, spkRates(j,:), x + dorsalEvents);  % actual
                if hasYhat; predicted{k,1,j} = interp1(fitdata.t, fitdata.yhat, x + dorsalEvents); end  % predicted
            end
            
            ventralEvents = data{i, 'ventral_touches'}{k};
            if ~isempty(ventralEvents)
                responses{k,2,j} = interp1(timeStamps, spkRates(j,:), x + ventralEvents);  % actual
                if hasYhat; predicted{k,2,j} = interp1(fitdata.t, fitdata.yhat, x + ventralEvents); end  % predicted
            end
        end
    end
    data{i, 'responses'} = {responses};
    data{i, 'predicted'} = {predicted};
end
fprintf('all done!\n')
% unitsPerSession = cellfun(@(x) size(x,3), data.responses);




%% plot that ish















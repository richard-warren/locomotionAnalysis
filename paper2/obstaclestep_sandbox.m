%% assess tuning steps over obstacle...


%% collect data


% settings
minTouchFrames = 4;  % only include touches in contact with obstacle for this number of consecutive frames (this is actually an approximation only... i use a simple medfilt for this, which is close to but not the same as duration thresholding)
x = linspace(-.4,.8,200);  % x axis for PSTHs

sessions = getEphysSessions();
n = length(sessions);

cellInits = repmat({{}, {}, {}, {}}, n, 1);
data = table(cellInits, cellInits, cell(n,1), cell(n,1), cell(n,1), 'RowNames', sessions, 'VariableNames', ...
    {'dorsal_touches', 'ventral_touches', 'responses', 'predicted', 'unit_ids'});

fprintf('preparing session: ')
for i = 1:length(sessions)
    fprintf('%i ', i)

    % load session touch times
    load(fullfile(getenv('OBSDATADIR'), 'sessions', sessions{i}, 'runAnalyzed.mat'), ...
        'frameTimeStamps', 'touchesPerPaw', 'touchClassNames', 'obsOnTimes', 'obsOffTimes');
    dorsalInds = find(contains(touchClassNames, 'orsal'));
    ventralInds = find(contains(touchClassNames, 'entral'));
    dorsalTouches = ismember(touchesPerPaw, dorsalInds);
    ventralTouches = ismember(touchesPerPaw, ventralInds);
    
    % debounce
    dorsalTouches = logical(medfilt1(double(dorsalTouches), minTouchFrames*2-1));
    ventralTouches = logical(medfilt1(double(ventralTouches), minTouchFrames*2-1));
    
    % for each contact type keep only first contact within each trial
    data{i, 'dorsal_touches'} = repmat({nan(length(obsOnTimes), 1)}, 1, 4);
    data{i, 'ventral_touches'} = repmat({nan(length(obsOnTimes), 1)}, 1, 4);
    for j = 1:length(obsOnTimes)
        bins = frameTimeStamps>obsOnTimes(j) & frameTimeStamps<obsOffTimes(j);
        dsub = dorsalTouches(bins,:);
        vsub = ventralTouches(bins,:);
        tsub = frameTimeStamps(bins);
        for k = 1:4
            dtouches = tsub(find(diff(dsub(:,k))==1, 1, 'first')+1);
            vtouches = tsub(find(diff(vsub(:,k))==1, 1, 'first')+1);
            if ~isempty(dtouches); data{i, 'dorsal_touches'}{k}(j) = dtouches; end
            if ~isempty(vtouches); data{i, 'ventral_touches'}{k}(j) = vtouches; end
        end
    end
    
    for k = 1:4  % get rid of nans
        data{i, 'ventral_touches'}{k}(isnan(data{i, 'ventral_touches'}{k})) = [];
        data{i, 'dorsal_touches'}{k}(isnan(data{i, 'dorsal_touches'}{k})) = [];
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
    data{i, 'unit_ids'} = {unit_ids};
end
unitsPerSession = cellfun(@(x) size(x,3), data.responses);
nunits = sum(unitsPerSession);
tic; save(fullfile('Y:\loco\obstacleData\data_transfer\to_remote', 'pawcontact.mat'), 'data', 'x'); toc  % save to remote storage
fprintf('all done!\n')




%% plot that ish

% settings
ncols = 24;  % divisible by 8 ideally
nrows = 14;
xlims = [-.2 .5];
colors = lines(8);  % colors for each contact type


close all
figure('name', 'paw_contact_responses', 'color', 'white', 'position', [1.00 1.00 2560.00 1363.00]); hold on
labels = {'LHd', 'LHv', 'LFd', 'LFv', 'RFd', 'RFv', 'RHd', 'RHv'};
offset = 0;

tic
for i = 1:length(sessions)  % for sessions
    disp(i/length(sessions))

    responses = data{i, 'responses'}{1};  % (paw) X (dorsal ventral) X (neuron)
    predicted = data{i, 'predicted'}{1};
    
    for j = 1:size(responses, 3)  % for units
        unitStartPlot = (sum(unitsPerSession(1:i-1))+j-1) * 8 + 1;  % start plot for this neuron
        
        for k = 1:4  % for paw
            for m = 1:2  % for dorsal, ventral
                plotInd = unitStartPlot + (k-1)*2 + m - 1 - offset;
                
                % make new figure if necessary
                if plotInd > ncols*nrows
                    offset = offset + ncols*nrows;
                    plotInd = plotInd - ncols*nrows;
                    saveas(gcf, ['E:\lab_files\paper2\plots\paw_contacts\psths' num2str(offset/(ncols*nrows)) '.png'])
                    figure('name', 'paw_contact_responses', 'color', 'white', 'position', [1.00 1.00 2560.00 1363.00]); hold on
                    disp('making new fig...')
                end
                
                subplot(nrows, ncols, plotInd); hold on  % dorsal
                cind = (k-1)*2 + m;
                
                if ~isempty(responses{k,m,j})
                    n = sum(~isnan(responses{k,m,j}(:,1)));
                    
                    % predicted
                    if ~isempty(predicted{k,m,j})
                        mean = nanmean(predicted{k,m,j}, 1);
                        stdev = nanstd(predicted{k,m,j}, 1);
                        patch([x fliplr(x)], [(-stdev+mean) fliplr(stdev+mean)], [.2 .2 .2], ...
                            'FaceAlpha', .25, 'EdgeColor', 'none')  % shaded error bars
                        plot(x, mean, 'color', [.2 .2 .2], 'lineWidth', 2)  % mean
                    end
                    
                    % response
                    mean = nanmean(responses{k,m,j}, 1);
                    if n>1
                        stdev = nanstd(responses{k,m,j}, 1);
                        patch([x fliplr(x)], [(-stdev+mean) fliplr(stdev+mean)], colors(cind,:), ...
                            'FaceAlpha', .25, 'EdgeColor', 'none')  % shaded error bars
                    end
                    plot(x, mean, 'color', colors(cind,:), 'lineWidth', 2)  % mean
                    
                    % vertical line at zero
                    plot([0 0], ylim, 'color', [.2 .2 .2])
                else
                    n = 0;
                end
                
                title(sprintf('%s %i', ...
                    sessions{i}, data{i, 'unit_ids'}{1}(j)), ...
                    'interpreter', 'none', 'FontSize', 8, 'FontWeight', 'normal')
                yticks = get(gca, 'ytick');
                set(gca, 'XLim', xlims, 'XTick', [], 'YTick', yticks([1,end]))
                
                xlabel(sprintf('%s (n=%i)', labels{cind}, n))
            end
        end
    end
end
set(gca, 'xcolor', get(gca, 'ycolor'))  % time axis visible only for final subplot
toc








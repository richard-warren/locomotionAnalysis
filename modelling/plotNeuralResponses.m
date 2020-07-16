function plotNeuralResponses(session, varargin)

% loads modelling\responses.mat and plots responses for all predictors //
% saves in obstacleData\figures\modelling\responses\

% settings
colors = lines(3);
s.eventColor = colors(1,:);
s.epochColor = colors(2,:);
s.contColor = colors(3,:);
s.showImportance = true;  % whether to show mutual information as text in figure
s.visible = true;  % whether figure is visible


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end % reassign settings passed in varargin
folder = fullfile(getenv('OBSDATADIR'), 'sessions', session);
load(fullfile(folder, 'modelling', 'predictors.mat'), 'predictors');
load(fullfile(folder, 'modelling', 'responses.mat'), 'responses');
if s.showImportance; load(fullfile(folder, 'modelling', 'importance.mat'), 'importance'); end
load(fullfile(folder, 'neuralData.mat'), 'unit_ids', 'spkRates', 'timeStamps');
fprintf('%s: plotting neural responses...\n', session)


% for each cell
for i = 1:length(unit_ids)
    
    
    % check that 'predictors' and 'responses' have same rows
    if ~isequal(predictors.Properties.RowNames, responses.Properties.RowNames)
        disp('WARNING! predictors and responses have different row names!')
        return;
    end
    
    % initialize figure
    fig = figure('color', 'white', 'name', sprintf('%s_cell%i', session, unit_ids(i)), ...
        'position', [185.00 58.00 1583.00 915.00], 'Visible', s.visible);
    rows = ceil(sqrt(height(responses)));  % same num row and cols
    yMax = prctile(spkRates(i,:), 99);

    
    for j = 1:height(responses)
        subplot(rows, rows, j); hold on
        name = responses.Properties.RowNames{j};
        xlabel(name, 'Interpreter', 'none')
        
        if responses.include(j)
            xLims = responses.xLims(j,:);

            if responses.type(j)=='event'
                plot([0 0], [0 yMax], 'color', [0 0 0 .4])
                response = responses.response{j}(:,:,i);
                x = linspace(xLims(1), xLims(2), size(response,2));
                
                respMean = nanmean(response,1);
                respStd = nanstd(response,1);
                plot(x, respMean, 'LineWidth', 3, 'color', s.eventColor)
                plot(x, respMean + [respStd; -respStd], 'LineWidth', 1, 'color', [s.eventColor .4])

            elseif responses.type(j)=='epoch'
                plot([0 0; 1 1]', [0 yMax; 0 yMax]', 'color', [0 0 0 .4])
                response = responses.response{j}(:,:,i);
                x = linspace(xLims(1), xLims(2), size(response,2));
                
                respMean = nanmean(response,1);
                respStd = nanstd(response,1);
                plot(x, respMean, 'LineWidth', 3, 'color', s.epochColor)
                plot(x, respMean + [respStd; -respStd], 'LineWidth', 1, 'color', [s.epochColor .4])

            elseif responses.type(j)=='continuous'
                response = responses.response{j}(:,i)';
                x = linspace(xLims(1), xLims(2), size(response,2));
                
                density = responses.density{j};
                spkRate = interp1(timeStamps, spkRates(i,:), predictors.t{j});
                density = density * (yMax/max(density));
                fill([xLims(1) x xLims(2) xLims(1)]', ...
                    [0 density 0 0]', [0 0 0], ...
                    'EdgeColor', [1 1 1]*.6, 'FaceAlpha', .1)
                inds = randperm(length(spkRate), 1000);
                scatter(predictors.data{j}(inds), spkRate(inds), 5, s.contColor, 'filled', 'markerfacealpha', .25)
                plot(x, response, 'LineWidth', 3, 'color', s.contColor)
            end
            
            set(gca, 'xlim', xLims, 'ylim', [0 yMax])

            % plot mutual information
            if s.showImportance
                mi = importance{name,'mi'}(i);
                text(xLims(2), yMax, sprintf('%.2f', mi), ...
                    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
            end
        end
    end
    
    % save figure
    pause(.1)
    saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures', 'modelling', 'responses', ...
        sprintf('%s cell%i responses.png', session, unit_ids(i))));
    if ~s.visible; close(fig); end
end


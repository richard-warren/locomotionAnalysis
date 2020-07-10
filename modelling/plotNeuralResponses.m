function plotNeuralResponses(session, varargin)

% loads modelling\responses.mat and plots responses for all predictors //
% saves in obstacleData\figures\modelling\responses\

% settings
colors = lines(3);
s.eventColor = colors(1,:);
s.epochColor = colors(2,:);
s.contColor = colors(3,:);


% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end % reassign settings passed in varargin
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'modelling', 'predictors.mat'), 'predictors');
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'modelling', 'responses.mat'), 'responses');
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'neuralData.mat'), ...
    'unit_ids', 'spkRates', 'timeStamps');



% for each cell
for i = 1:length(unit_ids)
    fprintf('%s, cell %i: plotting neural responses...\n', session, unit_ids(i))
    
    % check that 'predictors' and 'responses' have same rows
    if ~isequal(predictors.Properties.RowNames, responses(i).responses.Properties.RowNames)
        disp('WARNING! predictors and responses have different row names!')
        break;
    end
    
    % initialize figure
    cellResponses = responses(i).responses;
    figure('color', 'white', 'name', sprintf('%s_cell%i', session, unit_ids(i)), ...
        'position', [185.00 58.00 1583.00 915.00]);
    rows = ceil(sqrt(height(cellResponses)));  % same num row and cols
    yMax = prctile(spkRates(i,:), 99);

    
    for j = 1:height(cellResponses)
        
        response = cellResponses.response{j};
        name = cellResponses.Properties.RowNames{j};
        xLims = cellResponses.xLims(j,:);
        x = linspace(xLims(1), xLims(2), size(response,2));
        
        if cellResponses.type(j)=='event'
            subplot(rows, rows, j); hold on
            plot([0 0], [0 yMax], 'color', [0 0 0 .4])
            respMean = nanmean(response,1);
            respStd = nanstd(response,1);
            plot(x, respMean, 'LineWidth', 3, 'color', s.eventColor)
            plot(x, respMean + [respStd; -respStd], 'LineWidth', 1, 'color', [s.eventColor .4])
            xlabel(name, 'Interpreter', 'none')
            set(gca, 'xlim', xLims, 'ylim', [0 yMax])
            
        elseif cellResponses.type(j)=='epoch'
            subplot(rows, rows, j); hold on
            plot([0 0; 1 1]', [0 yMax; 0 yMax]', 'color', [0 0 0 .4])
            respMean = nanmean(response,1);
            respStd = nanstd(response,1);
            plot(x, respMean, 'LineWidth', 3, 'color', s.epochColor)
            plot(x, respMean + [respStd; -respStd], 'LineWidth', 1, 'color', [s.epochColor .4])
            xlabel(name, 'Interpreter', 'none')
            set(gca, 'xlim', xLims, 'ylim', [0 yMax])
        
        elseif cellResponses.type(j)=='continuous'
            subplot(rows, rows, j); hold on
            density = cellResponses.density{j};
            spkRate = interp1(timeStamps, spkRates(i,:), predictors.t{j});
            
            density = density * (yMax/max(density));
            fill([xLims(1) x xLims(2) xLims(1)]', ...
                [0 density 0 0]', [0 0 0], ...
                'EdgeColor', [1 1 1]*.6, 'FaceAlpha', .1)
            inds = randperm(length(spkRate), 1000);
            scatter(predictors.data{j}(inds), spkRate(inds), 5, s.contColor, 'filled', 'markerfacealpha', .25)

            plot(x, response, 'LineWidth', 3, 'color', s.contColor)
            xlabel(name, 'Interpreter', 'none')
            if xLims(1)~=xLims(2)  % somewhat of a hack to avoid constant predictors, which occur when feature is not successfully tracked...
                set(gca, 'xlim', xLims, 'ylim', [0 yMax])
            end
        end
    end
    
    % save figure
    pause(.1)
    saveas(gcf, fullfile(getenv('OBSDATADIR'), 'figures', 'modelling', 'responses', ...
        sprintf('%s cell%i responses.png', session, unit_ids(i))));
end
disp('all done!')


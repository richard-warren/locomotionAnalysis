function [models, fitdata] = fitEpochGlm(session, unit, varargin)
% fits models trained only during specific epochs (e.g. peri reward period)
% and test generalization outside of that epoch. for each epoch type,
% models are trained and only that epoch and with that epoch removed, then
% generalization is tested on the disjoint time periods // k folds
% validation is used within the epoch, and model trained on all time within
% epoch used to test generalization outside epoch // full models also
% trained // note that reward bases should not be used for this analysis,
% as the idea is to test whether neural relationship to BEHAVIOR is
% consitent across epochs


% settings
s.epochs = table({'reward_all'}, [0 4], ...  % epochs are defined by events, and the time windows (s) surroudning them
                 'RowNames', {'postreward'}, ...
                 'VariableNames', {'event', 'window'});
s.lambdas = logspace(-8, -1, 40);    % ridge regression coefficients
s.folds = 5;                         % cross-validation folds
s.parallel = false;                  % whether crossval analyses are parallelized
s.verbose = true;
s.save = true;
s.colors = [.4 .4 .4; lines(2)];     % firing rate // full model // model trained outside of epoch
s.outputFileName = fullfile(...
    getenv('SSD'), 'paper2', 'modelling', 'glms', 'epoch_glms', ...
    [session '_cell_' num2str(unit) '_glm.mat']);
s.plot = true;
s.closeFig = false;



% inits
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % parse name-value pairs
if s.verbose; fprintf('%s: fitting models for neuron %i... ', session, unit); end

% load design matrix
load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'designMatrices', 'epochs', [session '_designMatrix.mat']), ...
    'dmat', 't', 'reward_all')
load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'predictors', [session '_predictors.mat']), 'predictors')

% load neuron
neuralData = load(fullfile(getenv('SSD'), 'paper2', 'modelling', 'neuralData', [session '_neuralData.mat']));
spkRate = interp1(neuralData.timeStamps, neuralData.spkRates(neuralData.unit_ids==unit,:), t);
inds = find(~isnan(spkRate),1,'first') : find(~isnan(spkRate),1,'last');  % valid inds for neuron
spkRate = spkRate(inds);
t = t(inds);
dt = t(2) - t(1);
spkTimes = neuralData.spkTimes{neuralData.unit_ids==unit};
if isempty(spkTimes)
    fprintf('%s: WARNING! No spikes for neuron %i! Something is wrong. Aborting...\n', session, unit);
    return
end

% prep data for model
dmat = dmat(inds,:);
X = table2array(dmat);
y = histcounts(spkTimes, [t-dt/2 t(end)+dt/2])';


% define cross-validation folds
r = reward_all(reward_all>t(1) & reward_all<t(end))';  % epochs spanning beginning of one reward to the end of the next
r = [t(1) r t(end)];
fold_id = nan(size(t));
[~,~,f] = histcounts(randperm(length(r)-1), linspace(1,length(r)-1,s.folds+1));  % every ind assigned an int on [1,k]
for i = 1:s.folds
    for j = find(f==i)                     % for each reward epoch in fold
        fold_id(t>=r(j) & t<=r(j+1)) = i;  % assign each time point within epoch to fold i
    end
end


% initialize table
% (for each epoch, train a model in and out of the epoch (modelin,
% modelout) // for each of those models, save deviance when evaluated
% inside and outside the epoch (e.g. modelin_devin and modelin_devout) //
% finally, we compute deviance for full model both in and outside the epoch
epochNames = s.epochs.Properties.RowNames;
n = length(epochNames);
models = table(cell(n+1,1), cell(n+1,1), ...
    nan(n+1,1), nan(n+1,1), nan(n+1,1), nan(n+1,1), ...
    nan(n+1, 1), nan(n+1, 1), nan(n+1, length(t)), ...
    nan(n+1, length(t)), nan(n+1, length(t)), cell(n+1, 1), ...
    'RowNames', [{'full'}, s.epochs.Properties.RowNames], ...
    'VariableNames', {'modelin', 'modelout', 'modelin_devin', 'modelin_devout', ...
    'modelout_devin', 'modelout_devout', 'full_in', 'full_out', 'bins', 'yhatin', 'yhatout', 'eventTimes'});


% full
model_full = fitModel(X, y, s.lambdas);
models{'full', 'modelin'}{1} = model_full;
models{'full', 'modelin_devin'} = ...  % there is no 'in' vs 'out' for the full model, so just store deviance and yhat in 'in' column
    cvdeviance(X, y, model_full, 'holdout', true, 'bestLambdaOnly', true);
models{'full', 'yhatin'} = exp(model_full.fit_preval)' / dt;

for i = 1:height(s.epochs)
    % identify epoch time bins
    epochName = s.epochs.Properties.RowNames{i};
    eventName = s.epochs{i, 'event'}{1};
    events = predictors{eventName, 'data'}{1};
    models{epochName, 'eventTimes'} = {events};
    epochs = events + s.epochs{i, 'window'};         % start and end time of each epoch
    bins = any(t>=epochs(:,1) & t<=epochs(:,2), 1);  % time bins for epoch
    models{epochName, 'bins'} = bins;
    
    % deviance for full model in and out of epoch
    models{epochName, 'full_in'} = cvdeviance(X, y, model_full, ...
        'holdout', true, 'bestLambdaOnly', true, 'timeBins', bins);
    models{epochName, 'full_out'} = cvdeviance(X, y, model_full, ...
        'holdout', true, 'bestLambdaOnly', true, 'timeBins', ~bins);
    
    % train in, then out of epoch models
    for j = [{bins, 'in'}', {~bins, 'out'}']
        b = j{1};  % in or out of epoch bins
        c = j{2};  % 'in' or 'out' string
        
        % train model
        model = fitModel(X, y, s.lambdas, b);
        models{epochName, ['model' c]}{1} = model;

        % generate predictions FOR ALL t (covering periods in and out of training data)
        models{epochName, ['yhat' c]}(b) = exp(model.fit_preval) / dt;     % predictions in epoch (training)
        yhat = glmnetPredict(model.glmnet_fit, X(~b,:), [], 'response');   % predictions out of epoch (generalize)
        models{epochName, ['yhat' c]}(~b) = yhat(:, model.lambda_min_id) / dt;

        % compute deviance for data in and out of training set
        holdout = all(b==bins);
        models{epochName, ['model' c '_devin']} = ...
            cvdeviance(X, y, model, 'holdout', holdout, 'bestLambdaOnly', true, 'timeBins', bins);
        models{epochName, ['model' c '_devout']} = ...
            cvdeviance(X, y, model, 'holdout', ~holdout, 'bestLambdaOnly', true, 'timeBins', ~bins);
    end
end




if s.verbose; disp('all done!'); end


% save some objects for convenience...
fitdata.t = t;
fitdata.y = y;
fitdata.yRate = spkRate;
fitdata.session = session;
fitdata.unit = unit;

% save
settings = s;
if s.save; save(s.outputFileName, 'models', 'fitdata', 'settings'); end


if s.plot
    fig = figure('color', 'white', 'position', [2.00 722.00 1600 height(s.epochs)*250]); hold on
    eventToShow = {'reward_all', 'whiskerContact', 'lick'};
    eventColors = lines(length(eventToShow));
    yhatfull = models{'full', 'yhatin'};  % predictions for full model

    for i = 1:height(s.epochs)

        yhatout = models{epochName, 'yhatout'};  % predictions for model trained outside of epochs

        % sample trace
        subplot(height(s.epochs), 5, (i-1)*4+(1:3)); hold on
        epochName = s.epochs.Properties.RowNames{i};

        events = predictors{s.epochs.event{i}, 'data'}{1};
        events = events(events>t(1) & events<t(end));  % restrict to events during recording for this unit

        mid = events(round(length(events)/2)) + s.epochs.window(i,:);  % middle epoch
        xlims = mid + [-1 1]*range(s.epochs.window(i,:)) * 1.0;  % show twice the window size to the left and right
        events = events(events>xlims(1) & events<xlims(2));
        bins = t>xlims(1) & t<xlims(2);

        plot(t(bins), spkRate(bins), 'color', s.colors(1,:));  % real firing rate
        plot(t(bins), yhatfull(bins), 'color', s.colors(2,:));     % trained outside of epoch
        plot(t(bins), yhatout(bins), 'color', s.colors(3,:));     % trained outside of epoch

        legend('actual', 'train all', 'train out', 'autoupdate', 'off')

        % shade in epochs
        ylims = ylim;
        Xs = (events + [s.epochs.window(i,:) fliplr(s.epochs.window(i,:))])';  % every column has verticies of one rectangle
        Ys = repmat([ylims(1) ylims(1) ylims(2) ylims(2)]', 1, length(events));
        obj = patch(Xs, Ys, s.colors(3,:), 'EdgeColor', 'none', 'FaceAlpha', .1);
        uistack(obj, 'bottom')

        % add time bar
        plot(xlims(1) + [0 1], ylims(1) + [0 0], 'color', get(gca, 'XColor'), 'LineWidth', 2)
        text(xlims(1), ylims(1), '1 second', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')

        % add licks, rewards, other?
        for j = 1:length(eventToShow)
            events = predictors{eventToShow{j}, 'data'}{1};
            events = events(events>t(1) & events<t(end));  % restrict to events during recording for this unit
            plot([events events], ylims, 'color', eventColors(j,:))
        end
        
        set(gca, 'XLim', xlims, 'YLim', ylims, 'box', 'off', 'xcolor', 'none', 'TickDir', 'out')
        

        % psth
        xlims = s.epochs.window(i,:) + range(s.epochs.window(i,:))/2 * [-1 1];  % add half of window range as padding on either side
        x = linspace(xlims(1), xlims(2), 200);
        sigs = table([spkRate; yhatfull; yhatout], ...
                     {'actual', 'full', ['no ' epochName]}', ...
                     'VariableNames', {'signal', 'name'});  % pack all signals into a table
        events = predictors{s.epochs.event{i}, 'data'}{1};
        events = events(events>t(1) & events<t(end));  % restrict to events during recording for this unit
        subplot(height(s.epochs), 5, (i-1)*4+4); hold on

        for j = 1:3
            resp = interp1(t, sigs{j, 'signal'}, events + x);
            mn = nanmean(resp,1);
            stdev = nanstd(resp,1);
            patch([x fliplr(x)], [(-stdev+mn) fliplr(stdev+mn)], s.colors(j,:), ...
                    'FaceAlpha', .25, 'EdgeColor', 'none')           % shaded error bars
            plot(x, mn, 'color', s.colors(j,:), 'lineWidth', 2)  % mean
        end
        plot([0 0], ylim, 'color', [.2 .2 .2])
        title(s.epochs.event{i}, 'Interpreter', 'none')
        xlabel('time from event (s)')
        set(gca, 'TickDir', 'out')
        
        
        % deviance
        subplot(height(s.epochs), 5, (i-1)*4+5); hold on
        mat = [models{epochName, 'full_out'} models{epochName, 'full_in'}
               models{epochName, 'modelin_devout'} models{epochName, 'modelin_devin'}
               models{epochName, 'modelout_devout'} models{epochName, 'modelout_devin'}];  % (full x in x out) // (out x in)
        barFancy(repmat(mat,1,1,2), ...
            'levelNames', {{'train all', 'train in', 'train out'}, {'test out', 'test in'}}, ...
            'showScatter', false, 'lineThickness', .1, ...
            'colors', repelem([s.colors(2,:); .2 .2 .2; s.colors(3,:)], 2, 1))
    end
    
    if s.save; saveas(gcf, [s.outputFileName(1:end-4) '.png']); end
    if s.closeFig; close(fig); end
end



function fit = fitModel(X, y, lambdas, bins)
    % fit model with k fold cross validations
    if nargin<4; bins = true(1,size(X,1)); end
    if length(lambdas)==1; lambdas = [0 lambdas]; end  % a hack, because cvglmnet requires multiple lambdas
    
    % temporary hack to ensure all folds are represented
    % (should really re-assign temporally discontinuous chunks to different folds)
    if length(unique(fold_id(bins)))==s.folds
        ids = fold_id(bins);
    else
        ids = [];
    end
    
    if any(all(X==0,1)); disp('WARNING! Column with all zero predictors!!!'); end  % check for all-zero predictors
    
    options = struct('alpha', 0, 'lambda', lambdas, 'standardize', true);
    fit = cvglmnet(X(bins,:), y(bins), 'poisson', options, [], s.folds, ids, s.parallel, true);
    fit.fit_preval = fit.fit_preval(:,fit.lambda_min_id);  % only keep predictions for best lambda (save disk space)
end


end




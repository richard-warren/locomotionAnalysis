%% load cell feature importance
load(fullfile(getenv('OBSDATADIR'), 'matlabData', 'modelling', 'aggregates.mat'), 'aggregates', 'cellInfo');
mi = cat(2, aggregates.mi{:})';  % (predictor X cells) matrix of mutual information
aggregates.aggregate


%% plot n best cells per predictor

folder = 'Z:\loco\obstacleData\figures\modelling\cellExamples';
nBest = 5;
predictors = aggregates.Properties.RowNames;
slowPredictors = {'velocity', 'bodyAngle', 'whiskerAngle', 'buttHeight', 'satiation'};

for i = 57:length(predictors)
    [~, sortInds] = sort(mi(i,:), 'descend');
    
    % predictor specific settings
    if contains(predictors{i}, 'stride'); epochLims = [0 1]; else; epochLims = [-.25 1.25]; end
    if ismember(predictors{i}, slowPredictors)
        traceLims = [-15 2];
        kernelFall = .1;
    else
        traceLims = [-2 1];
        kernelFall = .02;
    end
    
    try
        for j = 1:nBest
            session = cellInfo{sortInds(j), 'session'}{1};
            unit = cellInfo{sortInds(j), 'unit'};
            plotPredictor(session, unit, predictors{i}, ... can i really put anything i 
                'epochLims', epochLims, 'traceLims', traceLims, 'kernelFall', kernelFall)

            name = sprintf('%s, session %s, unit %i', predictors{i}, session, unit);
            set(gcf, 'name', name)
            saveas(gcf, fullfile(folder, [name '.png']));
            pause(.1)
        end
        close all
    catch
        fprintf('%s: problem!\n', predictors{i})
    end
end
%% load cell feature importance
load(fullfile(getenv('SSD'), 'modelling', 'aggregates', 'aggregates.mat'), 'aggregates', 'cellInfo');
mi = cat(2, aggregates.mi{:})';  % (predictor X cells) matrix of mutual information

% % find location of each cell
ephysInfo = readtable(fullfile(getenv('OBSDATADIR'), 'spreadSheets', 'ephysInfo.xlsx'));
[~, inds] = ismember(unique(cellInfo.session, 'stable'), ephysInfo.session);
nucleus = ephysInfo.target(inds);

%% plot n best cells per predictor

folder = 'Z:\loco\obstacleData\figures\modelling\cellExamples';
nBest = 5;
predictors = aggregates.Properties.RowNames;
slowPredictors = {'velocity', 'bodyAngle', 'whiskerAngle', 'buttHeight', 'satiation'};
% nucleiColors = [lines(3); .2 .2 .2]*.8;
nucleiColors = repmat([.15 .15 .15], 4, 1);
nuclei = {'fastigial', 'interpositus', 'dentate', ''};



for i = 43:length(predictors)
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
    
    for j = 1:20
        try
            session = cellInfo{sortInds(j), 'session'}{1};
            unit = cellInfo{sortInds(j), 'unit'};
            
            sesNucleus = ephysInfo.target{strcmp(ephysInfo.session, session)};
            color = nucleiColors(strcmp(nuclei, sesNucleus),:);
            
            plotPredictorResponse(session, unit, predictors{i}, ...
                'epochLims', epochLims, 'traceLims', traceLims, 'kernelFall', kernelFall, ...
                'color', color)

            name = sprintf('%s, session %s, unit %i', predictors{i}, session, unit);
            set(gcf, 'name', name)
            title(sesNucleus)
            saveas(gcf, fullfile(folder, [name '.png']));
            pause(.1)
        catch
            fprintf('%s: problem!\n', predictors{i})
        end
    end
    close all
end

%% plot predictors in a nice order
close all
session = '200621_000';
xLims = [1058 1080];

% define variable 'groups'
limb = {'paw1LH_phase', 'paw1LH_x', 'paw1LH_x_vel', 'paw1LH_y', 'paw1LH_y_vel', 'paw1LH_z', 'paw1LH_z_vel', ...
    'paw2LF_phase', 'paw2LF_x', 'paw2LF_x_vel', 'paw2LF_y', 'paw2LF_y_vel', 'paw2LF_z', 'paw2LF_z_vel', ...
    'paw3RF_phase', 'paw3RF_x', 'paw3RF_x_vel', 'paw3RF_y', 'paw3RF_y_vel', 'paw3RF_z', 'paw3RF_z_vel', ...
    'paw1LH_phase', 'paw1LH_x', 'paw1LH_x_vel', 'paw1LH_y', 'paw1LH_y_vel', 'paw1LH_z', 'paw1LH_z_vel'};
gross = {'velocity', 'bodyAngle', 'bodyAngle_vel', 'buttHeight', 'buttHeight_vel'};
facial = {'ear', 'ear_vel', 'jaw', 'jaw_vel', 'lick', 'nose', 'nose_vel'};
whisker = {'whiskerAngle', 'whiskerAngle_vel', 'whiskerContact'};
reward = {'reward_normal', 'reward_omission', 'reward_surprise'};
visual = {'light'};
auditory = {'obstacle'};
ramps = {'obsToWisk'};
cutaneous = {'paw1LH_stance', 'paw2LF_stance', 'paw3RF_stance', 'paw4RH_stance'};
predictionError = {'paw1LH_contact_dorsal', 'paw1LH_contact_ventral', ...
    'paw2LF_contact_dorsal', 'paw2LF_contact_ventral', ...
    'paw3RF_contact_dorsal', 'paw3RF_contact_ventral', ...
    'paw4RH_contact_dorsal', 'paw4RH_contact_ventral'};

predictorList = [limb, gross, facial, whisker, reward, visual, auditory, ramps, cutaneous, predictionError];
plotNeuralPredictors(session, 'predictorList', predictorList, 'xLims', xLims)


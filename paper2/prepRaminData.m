%% prepare data needed by ramin to continue the paper

folder = 'E:\lab_files\paper2\ramin_data\';

%% unit metadata

unitInfo = getUnitInfo;
unitInfo = removevars(unitInfo, {'nucleus_target'});
writetable(unitInfo, fullfile(folder, 'unit_info.csv'));

%% predictors

% right now just manually copy/paste from E:\lab_files\paper2\modelling\predictors

%% neural_data

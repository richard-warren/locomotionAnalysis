% function getNeuralResponses(session)

% for each predictor in predictors.mat, save (trial X xaxis) matrix of
% responses to that predictor // responses are handled differently for
% events, epochs, and continuous vars // for events, we store simple PSTH
% firing rates // for epochs, each epoch (which varies in length) is
% stretched over a common x axis // for continuous vars, TBD!!!

% temp
session = '999999_999';

% settings


% initializations
cellData = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'cellData.csv'));
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'neuralData.mat'), ...
    'unit_ids', 'spkTimes', 'spkRates', 'timeStamps');
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'modelling', 'predictors.mat'), ...
    'contPredictors', 'epochPredictors', 'eventPredictors');


for i = 1:length(unit_ids)
    
    responses = 
    
    
    
end


% ------
% EVENTS
% ------


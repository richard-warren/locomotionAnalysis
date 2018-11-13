% function plotCellQualityMetrics(session, cell)

% temp
session = '181002_002';
cell = 1;

% initializations
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'neuralData.mat'));
spkTimes = unitTimes{cell};
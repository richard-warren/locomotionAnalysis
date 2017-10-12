
% plot baseline data

% user settings
dataDir = 'C:\Users\Rick\Google Drive\columbia\obstacleData\';
mouse = 'run5';

% initializations
load([dataDir 'sessions\171010_000\run.mat'])
load([dataDir 'sessions\171010_000\runAnalyzed.mat'])


%% plot stuff

close all; figure; plot(wheelTimes, wheelPositions); pimpFig
hold on; scatter(rewardTimes, 200*ones(1,length(rewardTimes)))
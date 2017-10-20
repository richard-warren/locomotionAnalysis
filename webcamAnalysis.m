
dataDir = 'C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\sessions\';
session = '171019_002';

% load data
load([dataDir session '\frameTimeStamps.mat'], 'timeStamps');

camSpikeClock = timeStamps;
camData = dlmread([dataDir session '\run.csv']) / 1000;
camSysClock = camData(:,1);
webCamSysClock = dlmread([dataDir session '\webCam.csv']) / 1000; % convert from ms to s

% remove discontinuities
timeSteps = cumsum([0; diff(webCamSysClock)<0]);
webCamSysClock = webCamSysClock + timeSteps;

timeSteps = cumsum([0; diff(camSysClock)<0]);
camSysClock = camSysClock + timeSteps;


sysToSpike = polyfit(camSysClock, camSpikeClock, 1);
webCamSpikeClock = webCamSysClock * sysToSpike(1) + sysToSpike(2);


%%
close all; figure;
scatter(camSysClock, camSpikeClock, 10, 'filled');
hold on
scatter(webCamSysClock, webCamTimeStamps, 25, 'filled');
pimpFig
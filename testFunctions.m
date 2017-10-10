%% load file

file = 'C:\Users\Rick\Google Drive\columbia\rick\run.mat';
addpath('spikeAnalysis', 'externalCode')
load(file);

% setup parameters
whEncoderSteps = 2880; % 720cpr * 4
wheelRad = 95.25; % mm
obEncoderSteps = 1000; % 250cpr * 4
obsRad = 96 / (2*pi); % radius of timing pulley driving belt of obstacles platform
targetFs = 1000;

%% test wheel decoder

[wheelPositions, wheelTimes] = rotaryDecoder(whEncodA.times, whEncodA.level,...
                                                     whEncodB.times, whEncodB.level,...
                                                     whEncoderSteps, wheelRad, targetFs);

figure;
plot(wheelTimes, wheelPositions, 'linewidth', 3);
xlabel('time (s)')
ylabel('position (m)')
pimpFig;

%% test obstacle decoder (from rotary encoder)
[obsPositions, obsTimes] = rotaryDecoder(obEncodA.times, obEncodA.level,...
                                                     obEncodB.times, obEncodB.level,...
                                                     obEncoderSteps, obsRad, targetFs);

figure;
plot(obsTimes, obsPositions, 'linewidth', 3);
xlabel('time (s)')
ylabel('position (m)')
pimpFig;


%% test motor decoder (commands sent to motor)

motorTimes = step.times;
motorPositions = motorDecoder(stepDir.level, stepDir.times, motorTimes, targetFs);

figure;
plot(motorTimes, motorPositions, 'linewidth', 3);
xlabel('time (s)')
ylabel('position (m)')
pimpFig;

%% compute normalized wheel position

rewardTimes = reward.times(diff(reward.values>2)==1);

wheelPositionsNorm = positionRewardNormalize(wheelPositions, wheelTimes, rewardTimes);

figure;
plot(wheelTimes, wheelPositionsNorm, 'linewidth', 3);
xlabel('time (s)')
ylabel('position (m)')
pimpFig;

%% plot wheel, motor, and obstacle positions
    
motorOnTimes = step.times(logical([1 ((diff(step.times')>.25)==1)])); % times at which platform movement becomes engaged (regardless of whether obstacle is also engaged)

motorPosNorm = motorPositions;
obsPosNorm = obsPositions;

% offset motor and obstacle positions to match wheel position
for i=1:length(motorOnTimes)
    
    futureMotorInds = motorTimes >= motorOnTimes(i);
    futureObsInds = obsTimes >= motorOnTimes(i);
    
    wheelPos = wheelPositionsNorm(find(wheelTimes>motorOnTimes(i), 1, 'first'));
    motorPos = motorPosNorm(find(motorTimes>motorOnTimes(i), 1, 'first'));
    obsPos = obsPosNorm(find(obsTimes>motorOnTimes(i), 1, 'first'));
    
    motorPosNorm(futureMotorInds) = motorPosNorm(futureMotorInds) - motorPos + wheelPos;
    obsPosNorm(futureObsInds) = obsPosNorm(futureObsInds) - obsPos + wheelPos;
       
end

% plot
figure;

plot(wheelTimes, wheelPositionsNorm, 'linewidth', 3); hold on
scatter(motorTimes, motorPosNorm, 50, 'filled');
scatter(obsTimes, obsPosNorm, 50, 'filled');

xlabel('time (s)')
ylabel('position (m)')
pimpFig;

%% test timestamp decoder

timeStamps = csvread('C:\Users\Rick\Google Drive\columbia\rick\timestampTest.csv');
times = timestampDecoder(timeStamps);

figure;
plot(times);
pimpFig;


%% test interpData

close all; figure
scatter(timesItp, positionsItp, 10, 'filled'); hold on
scatter(times, positions, 10, 'filled'); pimpFig


%% test getVelocity

load('C:\Users\Rick\Google Drive\columbia\rick\file drop\171009_003\run.mat');
whEncoderSteps = 2880; % 720cpr * 4
wheelRad = 95.25; % mm
windowSize = .2;
targetFs = 1000;
rewardTimes = reward.times(diff(reward.values>2)==1);

[wheelPositions, wheelTimes] = rotaryDecoder(whEncodA.times, whEncodA.level,...
                                                     whEncodB.times, whEncodB.level,...
                                                     whEncoderSteps, wheelRad, targetFs);

wheelPositionsNorm = positionRewardNormalize(wheelPositions, wheelTimes, rewardTimes);
%%
vel = getVelocity(wheelPositions, windowSize, targetFs);


% close all;
figure;
yyaxis left;
plot(wheelTimes, wheelPositionsNorm);
hold on

yyaxis right;
plot(wheelTimes, vel);

set(gca, 'xlim', [120 180])
pimpFig



























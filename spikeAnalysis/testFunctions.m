%% load file

file = 'C:\Users\LindseyBuckingham\Google Drive\columbia\analysis\run.mat';
load(file);

% setup parameters
whEncoderSteps = 2880; % 720cpr * 4
wheelRad = 95.25; % mm
obEncoderSteps = 1000; % 250cpr * 4
obsRad = 96 / (2*pi); % radius of timing pulley driving belt of obstacles platform

%% test wheel rotary decoder

[wheelPositions, wheelPositionTimes] = rotaryDecoder(whEncodA.times, whEncodA.level,...
                                                     whEncodB.times, whEncodB.level,...
                                                     whEncoderSteps, wheelRad);


close all;
figure;
plot(wheelPositionTimes, wheelPositions, 'linewidth', 3);
xlabel('time (s)')
ylabel('position (m)')
pimpFig;

%% test obstacle rotary decoder
[obsPositions, obsPositionTimes] = rotaryDecoder(obEncodA.times, obEncodA.level,...
                                                     obEncodB.times, obEncodB.level,...
                                                     obEncoderSteps, obsRad);


% close all;
figure;
plot(obsPositionTimes, obsPositions, 'linewidth', 3);
xlabel('time (s)')
ylabel('position (m)')
pimpFig;


%% test stepper decoder

motorPositions = motorDecoder(stepDir.level, stepDir.times, step.times);


% close all;
figure;
plot(step.times, motorPositions, 'linewidth', 3);
xlabel('time (s)')
ylabel('position (m)')
pimpFig;

%% test position reward normalization

posNorm = positionRewardNormalize(wheelPositions, wheelPositionTimes, reward.times);

close all;
figure;
plot(wheelPositionTimes, posNorm, 'linewidth', 3);
xlabel('time (s)')
ylabel('position (m)')
pimpFig;

%% plot wheel and motor position on top of one another

% offset motor position to match wheel position
obsTimes = obsOn.times;
obsTimes = obsTimes(logical(obsOn.level)); % only save inds where the obs is turning from off to on
motorPosNorm = motorPositions;

for i=1:length(obsTimes)
    
    inds = step.times>obsTimes(i);
    wheelPos = posNorm(find(wheelPositionTimes>obsTimes(i), 1, 'first'));
    motorPos = motorPosNorm(find(step.times>obsTimes(i), 1, 'first'));
    motorPosNorm(inds) = motorPosNorm(inds) - motorPos + wheelPos;
   
    
end


close all;
figure;
plot(wheelPositionTimes, posNorm, 'linewidth', 3); hold on
scatter(step.times, motorPosNorm, 5);
xlabel('time (s)')
ylabel('position (m)')
pimpFig;






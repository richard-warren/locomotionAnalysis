%% load file

file = 'C:\Users\Rick\Google Drive\columbia\rick\run.mat';
load(file);

% setup parameters
whEncoderSteps = 2880; % 720cpr * 4
wheelRad = 95.25; % mm
obEncoderSteps = 1000; % 250cpr * 4
obsRad = 96 / (2*pi); % radius of timing pulley driving belt of obstacles platform

%% test wheel decoder

[wheelPositions, wheelTimes] = rotaryDecoder(whEncodA.times, whEncodA.level,...
                                                     whEncodB.times, whEncodB.level,...
                                                     whEncoderSteps, wheelRad);

figure;
plot(wheelTimes, wheelPositions, 'linewidth', 3);
xlabel('time (s)')
ylabel('position (m)')
pimpFig;

%% test obstacle decoder (from rotary encoder)
[obsPositions, obsTimes] = rotaryDecoder(obEncodA.times, obEncodA.level,...
                                                     obEncodB.times, obEncodB.level,...
                                                     obEncoderSteps, obsRad);

figure;
plot(obsTimes, obsPositions, 'linewidth', 3);
xlabel('time (s)')
ylabel('position (m)')
pimpFig;


%% test motor decoder (commands sent to motor)

motorTimes = step.times;
motorPositions = motorDecoder(stepDir.level, stepDir.times, motorTimes);

figure;
plot(motorTimes, motorPositions, 'linewidth', 3);
xlabel('time (s)')
ylabel('position (m)')
pimpFig;

%% compute normalized wheel position

rewardTimes = reward.times(diff(reward.values>2)==1);

wheelPosNorm = positionRewardNormalize(wheelPositions, wheelTimes, rewardTimes);

figure;
plot(wheelTimes, wheelPosNorm, 'linewidth', 3);
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
    
    wheelPos = wheelPosNorm(find(wheelTimes>motorOnTimes(i), 1, 'first'));
    motorPos = motorPosNorm(find(motorTimes>motorOnTimes(i), 1, 'first'));
    obsPos = obsPosNorm(find(obsTimes>motorOnTimes(i), 1, 'first'));
    
    motorPosNorm(futureMotorInds) = motorPosNorm(futureMotorInds) - motorPos + wheelPos;
    obsPosNorm(futureObsInds) = obsPosNorm(futureObsInds) - obsPos + wheelPos;
       
end

% plot
figure;

plot(wheelTimes, wheelPosNorm, 'linewidth', 3); hold on
scatter(motorTimes, motorPosNorm, 50, 'filled');
scatter(obsTimes, obsPosNorm, 50, 'filled');

xlabel('time (s)')
ylabel('position (m)')
pimpFig;






%% load file

file = 'C:\Users\Rick\Google Drive\columbia\analysis\run.mat';
load(file);


%% test rotary decoder

[wheelPositions, positionTimes] = rotaryDecoder(encoderA.times, encoderA.level, encoderB.times, encoderB.level);


close all;
figure;
plot(positionTimes, wheelPositions, 'linewidth', 3);
xlabel('time (s)')
ylabel('position (m)')
pimpFig;



%% test stepper decoder

motorPositions = motorDecoder(stepDir.level, stepDir.times, step.times);


close all;
figure;
plot(step.times, motorPositions, 'linewidth', 3);
xlabel('time (s)')
ylabel('position (m)')
pimpFig;

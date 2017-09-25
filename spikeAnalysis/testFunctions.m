%% load file

file = 'C:\Users\Rick\Google Drive\columbia\analysis\run.mat';
load(file);


%% test rotary decoder

[positions, positionTimes] = rotaryDecoder(encoderA.times, encoderA.level, encoderB.times, encoderB.level);



close all;
figure;
plot(positionTimes, positions, 'linewidth', 3);
xlabel('time (s)')
ylabel('position (m)')
pimpFig;



%%

load('C:\Users\Rick\Google Drive\columbia\obstacleData\sessions\171014_002\run.mat')

%%
[positions, times] = rotaryDecoder(obEncodA.times, obEncodA.level, obEncodB.times, obEncodB.level, 1000, 95.25, 1000);
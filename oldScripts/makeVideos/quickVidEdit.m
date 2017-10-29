
% crops file from startTime to endTime


% user settings
startTime = 0;
endTime = 10;
fs = 250;
file = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\testVideo\topTest.mp4';

% initializations
frames = (startTime*fs)+1 : endTime*fs;

vid = VideoReader(file);
vidWrite = VideoWriter([file(1:(find(file=='.',1,'last')-1)) 'Edited.mp4'], 'MPEG-4');
open(vidWrite);


%% write selected frames to video
for i=1:length(frames)
   
    frame = read(vid, frames(i));
    writeVideo(vidWrite, frame);
    disp(i/length(frames));
    
end

close(vidWrite);
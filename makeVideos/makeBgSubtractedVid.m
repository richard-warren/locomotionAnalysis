function makeBgSubtractedVid(session)


% initializations
vid = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4']);
vidWriter = VideoWriter([getenv('OBSDATADIR') 'sessions\' session '\runBotBgSub.mp4'], 'MPEG-4');
open(vidWriter)
frameNum = vid.NumberOfFrames;
bg = getBgImage(vid, 1000, 120, 2*10e-4, false);


w = waitbar(0, 'creating background subtracted video...');

for i = 1:frameNum
    
    frame = read(vid, i) - bg;
    writeVideo(vidWriter, frame);
    
    waitbar(i/frameNum)
end

close(vidWriter)
close(w)
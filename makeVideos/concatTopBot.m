function concatTopBot(session)

vidTop = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);
vidBot = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4']);
vidWriter = VideoWriter([getenv('OBSDATADIR') 'sessions\' session '\runTopBot.mp4'], 'MPEG-4');
open(vidWriter)


w = waitbar(0, 'concatenating top and bottom views...');

for i = 1:vidTop.NumberOfFrames
    
    frameTop = read(vidTop, i);
    frameBot = read(vidBot, i);
    frame = cat(1, frameTop, frameBot);
    writeVideo(vidWriter, frame);
    
    waitbar(i/vidTop.NumberOfFrames)
end

close(vidWriter)
close(w)
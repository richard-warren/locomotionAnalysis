
% settings
session = '180615_005';
frameRate = 25;
obsPrePost = [-.05 .03];


vid = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runWisk.mp4']);
load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'],...
            'obsPositions', 'obsTimes', 'obsPixPositions', 'frameTimeStamps', 'obsOnTimes', 'obsOffTimes', 'nosePos');
obsPositions = fixObsPositions(obsPositions, obsTimes, obsPixPositions, frameTimeStamps, obsOnTimes, obsOffTimes, nosePos(1));
obsPosInterp = interp1(obsTimes, obsPositions, frameTimeStamps);


vidWriter = VideoWriter([getenv('OBSDATADIR') 'editedVid\' sprintf('wiskEg%s', session)], 'MPEG-4');
set(vidWriter, 'FrameRate', frameRate);
open(vidWriter);

for i = 1:length(obsOnTimes)
    disp(i/length(obsOnTimes))
    
    startInd = find(frameTimeStamps>obsOnTimes(i) & obsPosInterp>obsPrePost(1), 1, 'first');
    stopInd = find(frameTimeStamps<obsOffTimes(i) & obsPosInterp<obsPrePost(2), 1, 'last');
    
    for j = startInd:stopInd
        frame = read(vid, j);
        writeVideo(vidWriter, frame);
    end
end
close(vidWriter)
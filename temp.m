sessions = {'180122_001', '180122_002', '180122_003', ...
            '180123_001', '180123_002', '180123_003', ...
            '180124_001', '180124_002', '180124_003', ...
            '180125_001', '180125_002', '180125_003'};

        
for i = 1:length(sessions)
    
    vidBot = VideoReader([getenv('OBSDATADIR') 'sessions/' sessions{i} '/runBot.mp4']);
    
    [noseX, noseY, medianFrame] = getNosePos(vidBot);
    
    figure;
    imshow(medianFrame); hold on
    line([1 vidBot.Width], [noseY noseY], 'color', 'red', 'linewidth', 2);
    line([noseX noseX], [1 vidBot.Height], 'color', 'red', 'linewidth', 2);
    pimpFig
    
    
end




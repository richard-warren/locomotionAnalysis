% tracking example frames for running and licking

session = '200709_000';
scatSize = 80;
smatScatSizes = [5 25];
contrast = [.1 .7];
trial = 2;
color = [0.9290 0.6940 0.1250];
pawColors = lines(4);
removeRows = [200 210];
runInd = 2875;
lickInd = 800;
featuresToShow = {'paw1LH_top', 'paw2LF_top', 'paw3RF_top', 'paw4RH_top', ...
                  'paw1LH_bot', 'paw2LF_bot', 'paw3RF_bot', 'paw4RH_bot', ...
                  'nose_top', 'nose_bot', 'tailMid_top', 'tailMid_bot', ...
                  'tailBase_top', 'tailBase_bot', 'jaw', 'ear', 'tongue'};

close all
for ind = [runInd lickInd]
    showSingleFrameTracking(session, 1, 'ind', ind, ...
        'contrastLims', contrast, 'addWiskCam', true, 'pawColors', pawColors, ...
        'otherColors', color, 'featuresToShow', featuresToShow, ...
        'removeRows', removeRows, 'mainSize', scatSize, 'trailingSizes', smatScatSizes);
    saveas(gcf, fullfile(getenv('OBSDATADIR'), 'data_transfer', sprintf('%stracking_frame_%i.png', session, ind)))
    
    % save face cam only
    vidWisk = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runWisk.mp4'));
    frameWisk = 255 - read(vidWisk, ind);
    imwrite(frameWisk, fullfile(getenv('OBSDATADIR'), 'data_transfer', sprintf('%sface_frame_%i.png', session, ind)))
end


%%
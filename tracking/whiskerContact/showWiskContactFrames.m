function showWiskContactFrames(session)

% initializations
vid = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runWisk.mp4'));
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'wiskContactFrames', 'wiskContactFramesConfidences')


allImgs = zeros(vid.Height, vid.Width, length(wiskContactFrames));
for i = 1:length(wiskContactFrames)
    if wiskContactFrames(i)>0
        frame = read(vid, wiskContactFrames(i));
        frame = insertText(frame, [0 0], num2str(i), ...
            'FontSize', 40, 'BoxColor', [0 0 0], 'TextColor', 'white', 'BoxOpacity', 1);
        allImgs(:,:,i) = squeeze(frame(:,:,1));
    end
end

figure; montage(allImgs, 'DisplayRange', [0 255]);
set(gcf, 'name', session, 'units', 'normalized', 'position', [0 0 1 1])
function showWiskHeights(session)

% settings
bins = 6;  % number of height bins
trialsPerBin = 3;  % number of trials per bin

% initializations
fprintf('%s: showing whisker contact frames binned by obstacle height...\n', session)
vid = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runWisk.mp4'));
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'wiskContactFrames', 'wiskContactFramesConfidences', 'obsHeightsVid')

[hgts, edges] = discretize(obsHeightsVid, bins);

frames = zeros(vid.Height*bins, vid.Width*trialsPerBin);

for i = 1:bins
    inds = find(hgts==i);
    inds = sort(inds(randsample(length(inds), trialsPerBin)));

    for j = 1:length(inds)
        
        frame = read(vid, wiskContactFrames(inds(j)));
        
        frame = insertText(frame, [0 0], sprintf('%.1f', mean(edges(i:i+1))), ...
            'FontSize', 40, 'BoxColor', [0 0 0], 'TextColor', 'white', 'BoxOpacity', 1);
        
        ys = (i-1)*vid.Height+1:i*vid.Height;
        xs = (j-1)*vid.Width+1:j*vid.Width;
        frames(ys, xs) = frame(:,:,1);
    end
end

figure;
imshow(frames, [0 255])


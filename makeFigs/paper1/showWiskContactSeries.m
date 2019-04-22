function showWiskContactSeries(session, trials, framesPerRow)
% show a time series of frames in a row, with the flast frame being the
% detected contact frame

% settings
% session = '190319_002';
% trials = 10:10:30;
% framesPerRow = 6;
contrastLims = [.05 .4]; % pixels at these proportional values are mapped to 0 and 255
scaling = .5;
botCrop = 60; % trim this much off the bottom of the images
verticalSpacing = 20; % number of pixels separating rows


% initializations
vid = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runWisk.mp4'));
imgDims = [(vid.Height-botCrop)*length(trials) +  (length(trials)-1)*verticalSpacing, ...
           vid.Width*framesPerRow];
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), 'wiskContactFrames');



img = uint8(ones(imgDims)*255);
for trial = 1:length(trials)
    frames = read(vid, [wiskContactFrames(trials(trial))-framesPerRow+1, wiskContactFrames(trials(trial))]);
    frames = squeeze(frames(:,:,1,:)); % remove color channel
    row = 255 - reshape(frames, vid.Height, []);
    row = row(1:end-botCrop,:);
    y = (trial-1)*(vid.Height-botCrop) + 1 + (trial-1)*verticalSpacing;
    img(y:(y+vid.Height-botCrop-1),:) = row;
end
img = imadjust(img, contrastLims, [0 1]);


file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'wiskContactImgs', ...
        'whiskerContactImgs.png');
fprintf('writing %s to disk...\n', file)
imwrite(img, file);


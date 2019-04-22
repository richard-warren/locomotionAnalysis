function showPawContactSeries(session, preContactFrames, totalFrames)

% shows paw series of frames surrounding paw contacts of different types to
% demonstrate the contact analysis

% settings
% session = '180703_000';
% preContactFrames = 2; % how many frames to show before the contact occurs
% totalFrames = 8;
touchLengthLims = [2 inf]; % only include touches with consecutive frames between these limits
touchTypes = {'fore_ventral', 'fore_dorsal', 'hind_ventral_low', 'hind_dorsal'};
obsCropping = 100; % frames are cropped around obstacle in top view
contrastLims = [.1 .7]; % pixels at these proportional values are mapped to 0 and 255

verticalSpacing = 10; % number of pixels separating rows

% initializations
vid = VideoReader(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runTop.mp4'));
load(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'runAnalyzed.mat'), ...
    'touches', 'touchClassNames', 'touchesPerPaw');
obsTracking = readtable(fullfile(getenv('OBSDATADIR'), 'sessions', session, 'trackedFeaturesRaw.csv')); % get raw tracking data
obsTracking = table2array(obsTracking(:,contains(obsTracking.Properties.VariableNames, 'obs_top')));
[~, touchInds] = ismember(touchTypes, touchClassNames);
imgDims = [length(touchTypes)*obsCropping + (length(touchTypes)-1)*verticalSpacing, totalFrames*obsCropping];
img = uint8(ones(imgDims)*255);


for i = 1:length(touchTypes)
    
    t = touches==touchInds(i) & any(touchesPerPaw,2)'; % touches for features i
    
    % find touches that are of correct duration
    tStarts = find(diff(t)==1)+1;
    tEnds = find(diff(t)==-1);
    durations = tEnds - tStarts;
    tStarts = tStarts(durations>=touchLengthLims(1) & durations<touchLengthLims(2));
    
%     tStart = tStarts(randperm(length(tStarts))); % select random trial
    tStart = tStarts(1); % take first trial
    
    % get frames
    rowInds = tStart-preContactFrames : tStart-preContactFrames+totalFrames-1;
    row = uint8(ones(obsCropping,totalFrames*obsCropping)*255);
    for j = 1:length(rowInds)
        x = round(obsTracking(rowInds(j),1)) + .5*[-obsCropping obsCropping];
        y = round(obsTracking(rowInds(j),2)) + .5*[-obsCropping obsCropping];
        frame = rgb2gray(read(vid, rowInds(j)));
        try
            frame = frame(y(1):y(2)-1, x(1):x(2)-1); % get subframe
            if t(rowInds(j)); frame(end-2:end,:) = 255; end % underline images where contact is detected
            row(:, (j-1)*obsCropping+1 : j*obsCropping) = frame;
        end % put subframe in row
    end
    
    % add row to overall image
    y = (i-1)*(obsCropping+verticalSpacing)+1;
    img(y:y+obsCropping-1,:) = row;
end
img = imadjust(img, contrastLims, [0 1]); % adjust contrast


file = fullfile(getenv('OBSDATADIR'), 'papers', 'paper1', 'figures', 'pawContactImgs', ...
        'pageContactImg.png');
fprintf('writing %s to disk...\n', file)
imwrite(img, file);


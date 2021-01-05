

% settings
session = '200113_000'; %'180715_004';
overlayNum = 5;
imNum = 6;
cmap = 'hsv';
colorBump = .25;
nrows = 3;
contrast = [.075 .3];
% edgeFade = 20;
mask = [202 208];  % cut out these cols to remove bright spot between top and bottom views
filename = 'Z:\loco\obstacleData\data_transfer\title_graphic.png';

% inits
colors = eval([cmap '(overlayNum)']);
% vid = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\run.mp4']);
% load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], 'obsOnTimes', 'obsPixPositions', 'frameTimeStamps');
wid = vid.Width;
hgt = vid.Height - diff(mask) - 1;
obsLocations = fliplr(linspace(.1,.9,imNum) * wid); % pixel locations of obs

% fadeMask = [linspace(.5,1,edgeFade), ones(1,wid-2*edgeFade), linspace(1,.5,edgeFade)];
% fadeMask = repmat(fadeMask, vid.Height, 1);

imgs = cell(1, nrows);
for r = 1:nrows
    disp(r)
    trials = randperm(length(obsOnTimes), overlayNum);
    allImages = nan(hgt, wid, overlayNum, imNum);

    % get all images for all trials and obsLocations
    for i = 1:length(trials)
        for j = 1:length(obsLocations)
            frameInd = find(frameTimeStamps>=obsOnTimes(trials(i)) & obsPixPositions'<=obsLocations(j), 1, 'first');
            frame = double(rgb2gray(read(vid, frameInd)));
            frame(mask(1):mask(2), :) = [];            
%             frame = round(frame .* fadeMask);
            allImages(:,:,i,j) = frame;
        end
    end

    % create montage
    imMontage = (nan(hgt, wid*imNum, 3));
    for i = 1:imNum    
        colorIms = (nan(hgt, wid, 3, overlayNum));

        for j = 1:overlayNum
            img = allImages(:,:,j,i);
            imgColored = cat(3, img*colors(j,1), img*colors(j,2), img*colors(j,3));
            colorIms(:,:,:,j) = uint8(imgColored);
        end

        imMontage(:, (i-1)*wid+1:(i*wid), :) = sum(colorIms, 4);
    end

    % find grayness
    coloredness = mean(abs(diff(imMontage,1,3)),3) ./ (sum(imMontage,3)); % map of how colored vs. gray every pixel is
    imMontageColored = repmat(coloredness,1,1,3) .* imMontage; % imMontage where the pixels with color are accentuated
    imMontageMixed = imMontage*(1-colorBump)*mean(coloredness(:)) + imMontageColored*colorBump;
    imMontageHighRes = uint16(imMontageMixed .* ((2^16-1) / max(imMontageMixed(:))));
    imgs{r} = imMontageHighRes;
end

img = cat(1, imgs{:});
imgAdjusted = imadjust(img, contrast);
close all; figure;
imshow(imgAdjusted)
imwrite(imgAdjusted, filename)




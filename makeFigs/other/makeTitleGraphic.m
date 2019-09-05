

% settings
session = '180715_004';
overlayNum = 5;
imNum = 10;
colormap = 'hsv';
colorBump = .25;

% initializations
colors = eval([colormap '(overlayNum)']);
vidTop = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);
vidBot = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runBot.mp4']);
load([getenv('OBSDATADIR') 'sessions\' session '\runAnalyzed.mat'], 'obsOnTimes', 'obsPixPositions', 'frameTimeStamps');
obsLocations = fliplr(linspace(.1,.9,imNum) * vidTop.Width); % pixel locations of obs

trials = randperm(length(obsOnTimes), overlayNum);

allImages = nan(vidTop.Height+vidBot.Height, vidTop.Width, overlayNum, imNum);

% get all images for all trials and obsLocations
for i = 1:length(trials)
    for j = 1:length(obsLocations)
        frameInd = find(frameTimeStamps>=obsOnTimes(trials(i)) & obsPixPositions'<=obsLocations(j), 1, 'first');
        frameTop = rgb2gray(read(vidTop, frameInd));
        frameBot = rgb2gray(read(vidBot, frameInd));
        allImages(:,:,i,j) = cat(1, frameTop, frameBot);
    end
end



% create montage
imMontage = (nan(vidTop.Height+vidBot.Height, vidTop.Width*imNum, 3));

for i = 1:imNum    
    colorIms = (nan(vidTop.Height+vidBot.Height, vidTop.Width, 3, overlayNum));
    
    for j = 1:overlayNum
        img = allImages(:,:,j,i);
        imgColored = cat(3, img*colors(j,1), img*colors(j,2), img*colors(j,3));
        colorIms(:,:,:,j) = uint8(imgColored);
    end
    
    imMontage(:, (i-1)*vidTop.Width+1:(i*vidTop.Width), :) = sum(colorIms, 4);
end


% find grayness

coloredness = mean(abs(diff(imMontage,1,3)),3) ./ (sum(imMontage,3)); % map of how colored vs. gray every pixel is
imMontageColored = repmat(coloredness,1,1,3) .* imMontage; % imMontage where the pixels with color are accentuated
imMontageMixed = imMontage*(1-colorBump)*mean(coloredness(:)) + imMontageColored*colorBump;



imMontageHighRes = uint16(imMontageMixed .* ((2^16-1) / max(imMontageMixed(:))));
imwrite(imMontageHighRes, [getenv('OBSDATADIR') 'figures\imMontage.png'])

figure; imshow(imMontageHighRes)






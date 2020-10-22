

session = '180122_003';
vid = VideoReader([getenv('OBSDATADIR') 'sessions\' session '\runTop.mp4']);

%% phase bins

% settings
swingInds = [11996 12015];
binNum = 5;
yLims = [50 120];
xLims = [180 360];
contrastLims = [.2 1]; % pixels at these proportional values are mapped to 0 and 255
scaling = .8;

% initializations
frameInds = round(linspace(swingInds(1), swingInds(2), binNum+1));
frameInds = fliplr(frameInds(1:end-1));

dims = [diff(yLims)+1 diff(xLims)+1];
imgs = nan(dims(1), dims(2), binNum);


for i = 1:length(frameInds)
    img = rgb2gray(read(vid, frameInds(i)));
    img = img(yLims(1):yLims(2), xLims(1):xLims(2));
    img = imadjust(img, contrastLims, [0 1]);
    imgs(:,:,i) = img;
    
    imwrite(imresize(img, scaling), [getenv('OBSDATADIR') 'figures\binImages\phase\img' num2str(i) '.png']);
end

% close all; figure; imshow(uint8(reshape(imgs, dims(1), dims(2)*binNum, 1)))

%% speed bins
frameInd = 11999;
binNum = 5;
yLims = [50 120];
xLims = [200 330];
dims = [diff(yLims)+1 diff(xLims)+1];
offsetsMinMax = [.1 .5];
imCopies = 4;
scaling = .8;
dimming = .2;
leftFade = .5;

% initializations
offsetPixels = round(linspace(offsetsMinMax(1), offsetsMinMax(2), binNum) * dims(2));
velImgDims = [dims(1), dims(2) + max(offsetPixels)*(imCopies-1)];
xMid = round(velImgDims(2)/2);
brightnesses = linspace(dimming, 1, imCopies);

fadeVector = linspace(0, 1, dims(2)*leftFade);
fadeVector = [fadeVector ones(1, dims(2)-length(fadeVector))];
fadeImg = repmat(fadeVector, dims(1), 1);

img = rgb2gray(read(vid, frameInd));
img = img(yLims(1):yLims(2), xLims(1):xLims(2));
img = imadjust(img, contrastLims, [0 1]);

for i = 1:binNum
    
    velImg = nan(velImgDims);
    
    imWid = dims(2) + offsetPixels(i)*(imCopies-1);
    
    for j = 1:imCopies
        
        leftInd = round(xMid - (.5*imWid) + (j-1)*offsetPixels(i));
%         if j>1; keyboard; end
        velImg(:,leftInd:leftInd+dims(2)-1) = max(cat(3, ...
            velImg(:,leftInd:leftInd+dims(2)-1), double(img*brightnesses(j)).*fadeImg), [], 3);
        
    end
    
    imwrite(imresize(uint8(velImg), scaling), [getenv('OBSDATADIR') 'figures\binImages\speed\img' num2str(i) '.png']);
    
end


%% phase and speed bins
binNum = 5;
yLims = [50 120];
xLims = [160 360];
dims = [diff(yLims)+1 diff(xLims)+1];
offsetsMinMax = [.1 .2];
imCopies = 4;
scaling = .8;
dimming = .2;
leftFade = .5;

% initializations
offsetPixels = round(linspace(offsetsMinMax(1), offsetsMinMax(2), binNum) * dims(2));
velImgDims = [dims(1), dims(2) + max(offsetPixels)*(imCopies-1)];
xMid = round(velImgDims(2)/2);
brightnesses = linspace(dimming, 1, imCopies);

fadeVector = linspace(0, 1, dims(2)*leftFade);
fadeVector = [fadeVector ones(1, dims(2)-length(fadeVector))];
fadeImg = repmat(fadeVector, dims(1), 1);



for i = 1:binNum
    
    img = rgb2gray(read(vid, frameInds(i)));
    img = img(yLims(1):yLims(2), xLims(1):xLims(2));
    img = imadjust(img, contrastLims, [0 1]);
    velImg = nan(velImgDims);
    
    imWid = dims(2) + offsetPixels(i)*(imCopies-1);
    
    for j = 1:imCopies
        
        leftInd = round(xMid - (.5*imWid) + (j-1)*offsetPixels(i));
%         if j>1; keyboard; end
        velImg(:,leftInd:leftInd+dims(2)-1) = max(cat(3, ...
            velImg(:,leftInd:leftInd+dims(2)-1), double(img*brightnesses(j)).*fadeImg), [], 3);
        
    end
    
    imwrite(imresize(uint8(velImg), scaling), [getenv('OBSDATADIR') 'figures\binImages\phaseSpeed\img' num2str(i) '.png']);
    
end










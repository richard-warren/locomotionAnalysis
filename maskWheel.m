

% settings
thresh = 20;
lowerCircPts = [40 160; 216 115; 375 148];
upperCircPts = [24 158; 208 107; 371 136];
convKernel = [-1 0 1]';
circMinX = 40;

% initializations
midPointDiff = lowerCircPts(3,2) - lowerCircPts(2,2);
vidFile = 'C:\Users\LindseyBuckingham\Google Drive\columbia\obstacleData\svm\testVideo\runTop.mp4';
vid = VideoReader(vidFile);
lowerMask = getWheelMask(lowerCircPts, [vid.Height vid.Width]);
upperMask = getWheelMask(upperCircPts, [vid.Height vid.Width]);
wheelMask = uint8(lowerMask & ~upperMask);
[r, c] = fitCircle(lowerCircPts);
xs = circMinX:vid.Width;

figure;
imPreview = image(read(vid,1), 'CDataMapping', 'scaled');
colormap gray
set(gca, 'CLim', [0 255]);
pimpFig

for i = 1:vid.NumberOfFrames
    
    % get frame
    frame = read(vid,i);
    filtered = imfilter(rgb2gray(frame), convKernel, 'corr');
    filtered = filtered .* wheelMask;
    threshed = filtered>thresh;
    points = repmat(threshed, 1, 1, 3);
    
    % get circle coordinates
    ys = nan(1,length(xs));
    
    for j = 1:length(xs)
        col = threshed(:,xs(j));
        y = find(col);
        if isempty(y)
            ys(j) = nan;
        else
            ys(j) = median(y);
%             points(ys(j), xs(j), :) = [255 0 0];
        end
    end
    
    % fit circle
    xs = xs(~isnan(ys));
    ys = ys(~isnan(ys));
    [circR, circC] = circfit([xs', ys'], ones(length(xs),1));
    
    
    
    % update preview
    set(imPreview, 'CData', threshed*255);
    if exist('circs', 'var'); delete(circs); end
    circs = viscircles([circX, circY], circR);
%     set(imPreview, 'CData', uint8(threshed)*255);
                                            
    pause(.005);


end
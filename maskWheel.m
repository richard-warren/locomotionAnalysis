

% settings
thresh = 60;
startSearchInd = 115;
circRoiPts = [23 162; 209 109; 386 145];

% initializations
midPointDiff = circRoiPts(3,2) - circRoiPts(2,2);
vidFile = 'C:\Users\rick\Google Drive\columbia\obstacleData\svm\testVideo\runTop.mp4';
vid = VideoReader(vidFile);

figure;
imPreview = image(read(vid,1), 'CDataMapping', 'scaled');
colormap gray
set(gca, 'CLim', [0 255]);


for i = 1:vid.NumberOfFrames
    
    % get frame
    frame = rgb2gray(read(vid,i));
    threshed = frame>thresh;
    
    % update first and last circRoiPts
    circRoiPts(1,:) = [find(threshed(end,:), 1, 'first') vid.Height];
    circRoiPts(3,:) = [vid.Width find(threshed(startSearchInd:end,end), 1, 'first')+startSearchInd-1];
    circRoiPts(2,2) = circRoiPts(3,2) - midPointDiff;
    wheelMask = getWheelMask(circRoiPts, size(frame));
    
    % update preview
%     set(imPreview, 'CData', frame .* wheelMask);
    set(imPreview, 'CData', threshed*255);
    
    pause(.005);


end
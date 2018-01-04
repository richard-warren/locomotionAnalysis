function potentialLocationsTop = getPotentialLocationsTopMarkers(vid, frameInds, thresh, showTracking)


% !!! need to document


% settings
circRoiPts = [42 164; 220 115; 377 143];
circRoiOffSet = 4;
minBlobArea = 10;


% initializations
circRoiPts = circRoiPts - repmat([0 circRoiOffSet], 3, 1);
sampleFrame = rgb2gray(read(vid,1));
totalFrames = vid.NumberOfFrames;
wheelMask = getWheelMask(circRoiPts, [vid.Height vid.Width]);
bg = getBgImage(vid, 1000, 120, 2*10e-4, false);



% prepare figure
if showTracking

    fig = figure('menubar', 'none', 'color', 'white'); colormap gray

    rawAxis = subaxis(2,1,1, 'spacing', 0, 'margin', 0);
    rawIm = image(sampleFrame, 'parent', rawAxis, 'CDataMapping', 'scaled'); hold on;
    set(gca, 'visible', 'off');
    scatterPtsRaw = scatter(rawAxis, 0, 0, 50, 'filled', 'red');
    
    threshAxis = subaxis(2,1,2, 'spacing', 0, 'margin', 0);
    treshIm = image(sampleFrame, 'parent', threshAxis, 'CDataMapping', 'scaled'); hold on
    set(gca, 'visible', 'off');
    scatterPtsThresh = scatter(threshAxis, 0, 0, 50, 'filled', 'red');
    
    set(gcf, 'position', [1950 10 vid.Width*2 vid.Height*4])
    
end



% track markers on frameInd at a time
potentialLocationsTop(frameInds(end)) = struct();

for i = frameInds
    
    disp(i/totalFrames)
    
    % get frame and subframes
    frame = rgb2gray(read(vid,i)) .* wheelMask;
    frame = frame - bg;
    frameThreshed = frame > thresh;
    
    
    % blob analysis
    blobInfo = regionprops(frameThreshed, 'Area', 'Centroid'); % get blobs
    [~, sortInds] = sort([blobInfo.Area], 'descend');
    blobInfo = blobInfo(sortInds); % sort from largest to smallest
    if length(blobInfo) > 4; blobInfo = blobInfo(1:4); end % only keep four largest blobs
    blobInfo = blobInfo([blobInfo.Area] > minBlobArea); % get rid of blobs that are too small
            
    
    % store data
    centroids = reshape([blobInfo.Centroid], 2, length(blobInfo))';
    potentialLocationsTop(i).x = centroids(:,1);
    potentialLocationsTop(i).z = centroids(:,2);
    
    
    if showTracking
                
        % update figure
        set(rawIm, 'CData', frame);
        set(treshIm, 'CData', frameThreshed);
        set(scatterPtsRaw, 'XData', potentialLocationsTop(i).x, 'YData', potentialLocationsTop(i).z);
        set(scatterPtsThresh, 'XData', potentialLocationsTop(i).x, 'YData', potentialLocationsTop(i).z);
        
        % pause to reflcet on the little things...
        pause(.02);
%         keyboard
    end
end

if showTracking; close fig; end


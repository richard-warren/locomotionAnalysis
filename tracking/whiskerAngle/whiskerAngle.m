%% some brief tests on the ouput of janelia whisker tracking algo

% temp
addpath(genpath('C:\Program Files\WhiskerTracking\matlab'))
file = 'C:\Users\rick\Desktop\wiskTest\test';

% settings
minWiskLength = 150;  % pixels
minScore = 0;
showResults = true;
angleLims = [-40 140];
follicle_x_zLim = 100;  % only include whiskers with x follicle value +-follicle_x_zLim zscores


% initializations
wisks = LoadWhiskers([file '.whiskers']);
wisks = rmfield(wisks, {'id', 'thick', 'scores'}); % remove fields i don't need to free up memory
measures = LoadMeasurements([file '.measurements']);
measures = rmfield(measures, {'wid', 'label', 'face_x', 'face_y', 'curvature'}); % remove fields i don't need to free up memory
frameNum = double(wisks(end).time+1);
if showResults
    vid = VideoReader([file '.mp4']);
    close all; figure();
    im = imshow(read(vid,1)); hold on
    set(gcf, 'position', [1987 98 988 821])
%     bg = getBgImage(vid, 100, 0, 0, false);
end

% find valid whiskers
isLongEnough = [measures.length] > minWiskLength;
isValidScore = [measures.score] > minScore;
isValidAngle = ([measures.angle] > angleLims(1)) & ([measures.angle] < angleLims(2));
% isValidX = abs(zscore([measures.follicle_x])) < follicle_x_zLim & ...
%            abs(zscore([measures.tip_y])) < follicle_x_zLim;
isValid = isLongEnough & isValidScore;


for i = 1:frameNum
    
    disp(i/frameNum)
    pause(.05)
    
    % get frame and frameBins
    frame = read(vid,i);
    wiskBins = ([wisks.time] == i-1) & isValid;
    
    % get wisk data
    x = {wisks(wiskBins).x}; x = cat(1, x{:});
    y = {wisks(wiskBins).y}; y = cat(1, y{:});
    
    % ensure x and y acceptable frame indices
    x = round(x); x(x<1) = 1; x(x>size(frame,2)) = size(frame,2);
    y = round(y); y(y<1) = 1; y(y>size(frame,2)) = size(frame,2);
    
    % add tracked wisks to frame
    inds = sub2ind(size(frame), y, x);
    frame(inds) = 255;
    
    % update preview
    set(im, 'CData', frame);
    pause(.001)
    
end
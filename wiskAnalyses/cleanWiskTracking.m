
% temp
session = '171231_002';

% initializations
file = ['C:\Users\rick\Google Drive\columbia\obstacleData\sessions\' session '\runWiskEdited.'];
wisks = LoadWhiskers([file 'whiskers']);
measures = LoadMeasurements([file 'measurements']);
vid = VideoReader([file 'mp4']);


%%
close all; figure;
im = imshow(read(vid,1));


for i = 1:vid.NumberOfFrames
    
    % get frame
    frame = read(vid,i);
    
    % get wisk data
    bins = [wisks.time] == i-1; % measurements are zero indexed...
    x = {wisks(bins).x}; x = cat(1, x{:});
    y = {wisks(bins).y}; y = cat(1, y{:});
    
    % make x and y acceptable frame indices
    x = round(x); x(x<1) = 1; x(x>size(frame,2)) = size(frame,2);
    y = round(y); y(y<1) = 1; y(y>size(frame,2)) = size(frame,2);
    
    % add tracked wisks to frame
    inds = sub2ind(size(frame), y, x);
    frame(inds) = 0;
    
    % update preview
    set(im, 'CData', frame);
    pause(.1)
    
end
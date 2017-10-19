function getObsLocation(session, obsPos)

% displays all frames in a given session at which the obstacle is at position obsPos
% obstacle positions are based on the rotary encoder attached to the stepper motor track, and are corrected for drift using fixObsPositions
% it displays a montage of one frame per trial in which the obstacle is detected at obsPos
% it also displays the average of these frames // if the obstacle tracking is reliable all the frames should show the obstacle at the same location
% this can be used to determine a 'reference point' obstacle position that canbe used in other functions, ie plotting speed relative to position at which the obstacle is in the center of the wheel
%
% inputs     session:  name of session folder to be analyzed
%            obsPos:   obstacle position (m) relative to the very start of the obstacle track

% load data
dataDir = 'C:\Users\Rick\Google Drive\columbia\obstacleData\sessions\';
vid = VideoReader([dataDir session '\runTop.mp4']);

load([dataDir session '\runAnalyzed.mat'], 'obsPositions', 'obsTimes', 'obsOnTimes', 'obsOffTimes');
load([dataDir session '\frameTimeStamps.mat'], 'timeStamps');
obsPositions = fixObsPositions(obsPositions, obsTimes, obsOnTimes); % correct for drift in obstacle position readings

% user settings
windows = 3; % number of separate windows to plot all frames in

 

% collect all frames
allFrames = nan(vid.Height, vid.Width, 1, length(obsOnTimes));

for i = 1:length(obsOnTimes)
%     keyboard
    
    trialTime = obsTimes(find(obsTimes>obsOnTimes(i) & obsPositions >=obsPos, 1, 'first'));
    frameNum = find(timeStamps>=trialTime, 1, 'first');
    
    if length(frameNum)==1
        frame = rgb2gray(read(vid, frameNum));
        allFrames(:,:,:,i) = double(frame)/255;
    end
end


% plot all frames
groupInds = floor(linspace(1,length(obsOnTimes)+1, windows+1));

for i=1:windows
    
    figure;
    inds = groupInds(i):groupInds(i+1)-1;
    montage(allFrames(:,:,:,inds))
    pimpFig;
    
end


% plot average frame
figure; imshow(nanmean(allFrames,4)); pimpFig


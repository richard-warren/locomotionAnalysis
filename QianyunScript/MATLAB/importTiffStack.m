function xyzCoords = importTiffStack(filePath, downsampleRate, thickness)

% This function takes a mask tiff stack file and reformats it into a
% 3-dimension logical matrix which MATLAB can easily manipulate and handle.

% INPUT: filepath - the file path to the tiff stack you want to import.
%        downsampleRate - to what percentage you want to downsample your tiff stack imported. The default setting is downsampling to 30%.

% OUTPUT: a n-by-m-by-z three-dimension logical matrix. The n and m are length and width of the tiff image. The z is the number of tiffs in your tiff stack.
%         In the logical matrix, 0 means this point belongs to the background in your mask tiff file, and 1 means this point belongs to the feature
%         you traced in the mask tiff file. 



% settings
if ~exist('downsampleRate', 'var'); downsampleRate = 0.3; end
if ~exist('thickness', 'var'); thickness = 50; end % thickness of the brain sections

% import the mask tiff stack file
tiff_info = imfinfo(filePath);
tiff = imread(filePath, 1);
tiff_binary = tiff>0;
clear tiff
tiff_downsample = imresize(tiff_binary, downsampleRate);
clear tiff_binary
tiff_stack = false(size(tiff_downsample, 1), size(tiff_downsample, 2), length(tiff_info)); % initialization for the tiff stack
tiff_stack(:, :, 1) = tiff_downsample;
clear tiff_downsample

for i = 1 : length(tiff_info) 
    temp = imread(filePath, i);
    temp_binary = temp>0; % tuen the image into binary 
    clear temp
    temp_downsample = imresize(temp_binary, downsampleRate); 
    tiff_stack(:, :, i) = temp_downsample; 
    clear temp_downsample
end


xyzCoords = [];

resolution = 2; % The resolution of the tiff stack is 2um/pixel.

inds = find(tiff_stack>0);
s = size(tiff_stack);
[rows, cols, z] = ind2sub(s, inds);
xyzCoords = [cols*resolution/downsampleRate, (z-1)*thickness, rows*resolution/downsampleRate]; % covert pixel location into xyz coordinates (unit in microns)

clear tiff_stack


end
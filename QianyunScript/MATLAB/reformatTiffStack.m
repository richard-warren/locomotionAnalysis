function tiff_stack = reformatTiffStack(filePath, downsampleRate)

% This function takes a mask tiff stack file and reformats it into a
% 3-dimension logical matrix which MATLAB can easily manipulate and handle.

% INPUT: filepath - the file path to the tiff stack you want to import.
%        downsampleRate - to what percentage you want to downsample your tiff stack imported. The default setting is downsampling to 30%.

% OUTPUT: a n-by-m-by-z three-dimension logical matrix. The n and m are length and width of the tiff image. The z is the number of tiffs in your tiff stack.
%         In the logical matrix, 0 means this point belongs to the background in your mask tiff file, and 1 means this point belongs to the feature
%         you traced in the mask tiff file. 



% settings
if ~exist('downsampleRate', 'var'); downsampleRate = 0.3; end

% import the mask tiff stack file
tiff_info = imfinfo(filePath);
tiff = imread(filePath, 1); % read in first image
tiff_binary = tiff>0; % convert into binary 
clear tiff
tiff_downsample = imresize(tiff_binary, downsampleRate); % downsample the binary tiff
tiff_stack = false(size(tiff_downsample, 1), size(tiff_downsample, 2), length(tiff_info)); % initialization for the tiff stack

for i = 2 : length(tiff_info) 
    temp = imread(filePath, i);
    temp_binary = temp>0;
    clear temp
    temp_downsample = imresize(temp_binary, downsampleRate); 
    tiff_stack(:, :, i) = temp_downsample;    
end


end
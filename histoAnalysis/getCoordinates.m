function xyzCoords = getCoordinates(tiff, thickness, downsampleRate)

% This function is to get the xyz coordinates from the 3-dimension logical
% matrix created from the mask tiff file.




if ~exist('thickness', 'var'); thickness = 50; end % thickness of the brain sections
if ~exist('downsampleRate', 'var'); downsampleRate = 0.3; end % take the image downsample rate into consideration
xyzCoords = [];

resolution = 2; % The resolution of the tiff stack is 2um/pixel.

inds = find(tiff>0);
s = size(tiff);
[rows, cols, z] = ind2sub(s, inds);
xyzCoords = [cols*resolution/downsampleRate, rows*resolution/downsampleRate, (z-1)*thickness]; % covert pixel location into xyz coordinates (unit in microns)


end
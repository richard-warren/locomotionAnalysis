function d2 = distanceConverter(d1, shrinkageRate, downsampleRate)

% This function is to convert real world distance (in microns) to histo
% tiff image distances (in pixels), considering the downsample rate and
% tissue shrinkage rate.

% settings
resolution = 2;  % um/pixel

if ~exist('shrinkageRate', 'var'); shrinkageRate = 1.0; end 
if ~exist('downsampleRate', 'var'); downsampleRate = 0.3; end 

% converting
d2 = (d1.*resolution)./(shrinkageRate*downsampleRate); 

end
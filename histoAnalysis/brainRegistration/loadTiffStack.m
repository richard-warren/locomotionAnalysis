function tiffStack = loadTiffStack(fileName, varargin)
% quickly load tiff stack for a single channel into an (imageNum X height X width matrix)

% code modified from https://www.mathworks.com/matlabcentral/answers/105739-how-to-show-tiff-stacks

% settings
s.scaling = 1;  % shrink images by setting < 1
s.dataType = 'uint8';  % set to 'logical' to binarize images

% initializations
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
warning('off', 'imageio:tiffmexutils:libtiffWarning') % Suppress all the tiff warnings
tiffObj  = Tiff(fileName);
img = readImg();
[hgt, wid] = size(img, 1:2);
imgNum = length(imfinfo(fileName));

% collect images
fprintf('loading %s with %.2f scaling, image number: 1 ', fileName, s.scaling)
tiffStack = cast(zeros(imgNum, hgt, wid), s.dataType);
tiffStack(1,:,:) = img(:,:,1);
for i = 2:imgNum
    tiffObj.nextDirectory()
    img = readImg();
    fprintf('%i ', i)
    tiffStack(i,:,:) = img(:,:,1);
end
fprintf('all done!\n')
warning('on', 'imageio:tiffmexutils:libtiffWarning')


function img = readImg()
    img = cast(tiffObj.read(), s.dataType);
    if s.scaling~=1; img = imresize(img, s.scaling); end
end


end
function overlay = overlayImgs(imgs, opts)

% given a series of grayscale images (height X width X imgNum) and a vector of
% colors (imgNum X 3), returns the images overlayed, with each image tinted
% in color and then overlayed

% settings
s.colors = 'hsv';
s.contrastLims = [.2 .6]; % this is applied to each image before creating overlay
s.cutoff = 100; % normalize final image s.t. max value is the s.cutoff percentile of the pre-normalized image
s.projection = 'max'; % set to 'max' or 'mean' - determines how the images are combined

% reassign settings contained in opts
if exist('opts', 'var'); for i = 1:2:length(opts); s.(opts{i}) = opts{i+1}; end; end

% initializations
if ischar(s.colors); s.colors = eval([s.colors '(size(imgs,3))']); end % set colorspace if color is specified as a string


% color each image
imgsColored = nan(size(imgs,1), size(imgs,2), 3, size(imgs,4));
for i = 1:size(imgs,3)
    img = imgs(:,:,i);
    img = imadjust(img, s.contrastLims, [0 1]);
    img = double(img);
    imgsColored(:,:,:,i) = cat(3, img*s.colors(i,1), img*s.colors(i,2), img*s.colors(i,3));
end

% overlay images
if strcmp(s.projection, 'mean')
    overlay = max(imgsColored,[],4); % overlay images
elseif strcmp(s.projection, 'max')
    overlay = mean(imgsColored, 4); % overlay images
end

% normalize and convert to uint8
cutoff = prctile(overlay(:), s.cutoff);
overlay = min(overlay * (255/cutoff), 255);
overlay = uint8(overlay);

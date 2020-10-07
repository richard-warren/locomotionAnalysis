function plotLabels3D(labels, varargin)


% settings
s.downSampling = 1;             % (int, >=1) downsample volumes to speed things up 
s.maxSmps = 5e6;                % if > s.maxSmps in volumes, automatically downsample to prevent crashing
s.surfArgs = {};                % args pass to trisurf function (for boundary and convhull methods)
s.method = 'boundary';          % 'boundary', 'convhull', or 'contours'
s.apGrid = 1:size(labels,1);
s.mlGrid = 1:size(labels,3);
s.dvGrid = 1:size(labels,2);

% these settings only apply when 'method' is 'contours'
s.contourAlpha = .4;            
s.nLines = 40;
s.contourSmoothing = 3;
s.contourOpening = 3;


% inits
regionLabels = unique(labels(:));
regionLabels = regionLabels(regionLabels~=0)';
s.colors = lines(length(unique(labels(:)))-1);
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
hold on; axis equal
s.surfArgs = [{'FaceAlpha', .4, 'EdgeColor', [.2 .2 .2], 'EdgeAlpha', .4}, s.surfArgs];  % add defaults

% check if there are too many samples
smpNum = sum(labels(:)>0);
if smpNum > s.maxSmps
    newDownSampling = ceil(smpNum / s.maxSmps);
    if newDownSampling > s.downSampling
        s.downSampling = newDownSampling;
        fprintf('WARNING! Downsamping by %ix because too many samples detected...\n', s.downSampling);
    end
end

% downsample
if s.downSampling ~= 1
    ds = s.downSampling;
    s.apGrid = s.apGrid(1:ds:end);
    s.mlGrid = s.mlGrid(1:ds:end);
    s.dvGrid = s.dvGrid(1:ds:end);
    labels   = labels(1:ds:end, 1:ds:end, 1:ds:end);
end

% inits for contours setting
if strcmp(s.method, 'contours')
    [X,Y,Z] = meshgrid(s.mlGrid, s.apGrid, s.dvGrid);
    mlslice = linspace(1, s.mlGrid(end), s.nLines);
    apslice = linspace(1, s.apGrid(end), s.nLines);
    dvslice = []; %linspace(1, s.dvGrid(end), s.nLines);  % uncomment for sagittal slices as well
end


for i = regionLabels
    volume = labels==i;
    [apInds, dvInds, mlInds] = ind2sub(size(labels), find(volume));
    ap = s.apGrid(apInds); ml = s.mlGrid(mlInds); dv = s.dvGrid(dvInds);
    
    switch s.method
        case 'boundary'
            k = boundary(ap', ml', dv');
            trisurf(k, ml, ap, dv, 'FaceColor', s.colors(i,:), s.surfArgs{:});
            
        case 'convhull'
            k = convhull(ap, ml, dv, 'Simplify', true);
            trisurf(k, ml, ap, dv, 'FaceColor', s.colors(i,:), s.surfArgs{:});
            
        case 'contours'
            vol = double(volume);
            vol = imclose(vol, strel('sphere', s.contourOpening));
            vol = imfill(vol, 'holes');
            vol = smooth3(vol, 'box', s.contourSmoothing);
            contours = contourslice(X, Y, Z, permute(vol, [1 3 2]), dvslice, apslice, mlslice, [.5 .5], 'cubic');
            for j = 1:length({contours.EdgeAlpha})
                contours(j).EdgeAlpha = s.contourAlpha;
                contours(j).EdgeColor = s.colors(i,:);
            end
    end 
end

xlabel('ml'); ylabel('ap'); zlabel('dv');
set(gca, 'ZDir', 'reverse');
view(-45, 30)

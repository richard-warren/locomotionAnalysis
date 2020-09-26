function plotLabels3D(labels, varargin)


% settings
s.apGrid = 1:size(labels,1);
s.mlGrid = 1:size(labels,3);
s.dvGrid = 1:size(labels,2);
s.surfArgs = {};                % args pass to trisurf function


% inits
regionLabels = unique(labels(:));
s.colors = lines(length(unique(labels(:)))-1);
if exist('varargin', 'var'); for i = 1:2:length(varargin); s.(varargin{i}) = varargin{i+1}; end; end  % reassign settings passed in varargin
if ~iscell(labels); labels = {labels}; end
regionLabels = regionLabels(regionLabels~=0)';
hold on; axis equal


for i = 1:length(labels)
    for j = regionLabels
        [apInds, dvInds, mlInds] = ind2sub(size(labels{i}), find(labels{i}==j));
        ap = s.apGrid(apInds); ml = s.mlGrid(mlInds); dv = s.dvGrid(dvInds);
        k = boundary(ap', ml', dv');
%         k = convhull(ap, ml, dv, 'Simplify', true);
        trisurf(k, ml, ap, dv, ...
            'FaceColor', s.colors(j,:), 'FaceAlpha', .4, 'EdgeColor', [.2 .2 .2], 'EdgeAlpha', .4, ...
            s.surfArgs{:});
    end
end

xlabel('ml'); ylabel('ap'); zlabel('dv');
set(gca, 'ZDir', 'reverse');
view(-45, 30)
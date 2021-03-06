% settings
paws = [2 3];
animate = false;
trial = randperm(length(data),1);
azimuth = -45;
elevation = 20;
projectionDarkening = .6;
obsThickness = .003;
obsHeight = .009;
limsX = [-.1 .1];
% limsX = [min(data(trial).locations(:,1,i)) max(data(trial).locations(:,1,i))];
limsY = [-0.0381 0.0381];
limsZ = [-.01 .04];




obsX = [0 obsThickness]-.5*obsThickness;
obsY = limsY;
obsZ = [0 obsHeight];


close all; figure;
cmap = hsv(4);
xz = cell(1,length(paws));
xy = cell(1,length(paws));
xyz = cell(1,length(paws));

for i = paws
    
    % xz
    xz{i} = plot3(data(trial).locations(:,1,i), ones(1,size(data(trial).locations,1))*limsY(2), data(trial).locations(:,3,i), ...
        'color', cmap(i,:)*(1-projectionDarkening), 'linewidth', 3); hold on
    
    % xy
    xy{i} = plot3(data(trial).locations(:,1,i), data(trial).locations(:,2,i), zeros(1,size(data(trial).locations,1)), ...
        'color', cmap(i,:)*(1-projectionDarkening), 'linewidth', 3); hold on
    
    % xyz
    xyz{i} = plot3(data(trial).locations(:,1,i), data(trial).locations(:,2,i), data(trial).locations(:,3,i), ...
        'color', cmap(i,:), 'linewidth', 3); hold on
end

% add obs
% patch(zeros(1,4), [limsY(1) limsY(1) limsY(2) limsY(2)], [obsLimsZ(1) obsLimsZ(2) obsLimsZ(2) obsLimsZ(1)], [0 0 0])
% specify corners of shape
vertices = [obsX(1) obsY(1) obsZ(1)
            obsX(1) obsY(2) obsZ(1)
            obsX(1) obsY(2) obsZ(2)
            obsX(1) obsY(1) obsZ(2)
            obsX(2) obsY(1) obsZ(1)
            obsX(2) obsY(2) obsZ(1)
            obsX(2) obsY(2) obsZ(2)
            obsX(2) obsY(1) obsZ(2)];
% specify which corners to connect for each of 6 faces (6 sides of obs)
faces = [1 2 3 4
         5 6 7 8
         1 2 6 5
         4 3 7 8
         2 6 7 3
         1 5 8 4];

patch('Vertices', vertices, 'Faces', faces);

% pimp fig
daspect([1 1 1]); pimpFig
set(gca, 'visible', 'off', 'XLim', limsX, 'YLim', limsY, 'ZLim', limsZ, 'view', [azimuth elevation])
line([min(data(trial).locations(:,1,i)) max(data(trial).locations(:,1,i))], repmat(limsY(2), 1,2), [0 0], ...
    'color', [0 0 0], 'linewidth', 2) % x1
line([min(data(trial).locations(:,1,i)) max(data(trial).locations(:,1,i))], repmat(limsY(1), 1,2), [0 0], ...
    'color', [0 0 0], 'linewidth', 2) % x2
blackenFig


% animate
if animate
    for j = 1:size(data(trial).locations)
        for i = paws

            set(xz{i}, 'XData', data(trial).locations(1:j,1,i), ...
                       'YData', ones(1,j)*limsY(2), ...
                       'ZData', data(trial).locations(1:j,3,i));
            set(xy{i}, 'XData', data(trial).locations(1:j,1,i), ...
                       'YData', data(trial).locations(1:j,2,i), ...
                       'ZData', zeros(1,j));
            set(xyz{i},'XData', data(trial).locations(1:j,1,i), ...
                       'YData', data(trial).locations(1:j,2,i), ...
                       'ZData', data(trial).locations(1:j,3,i));
            pause(.02);
        end
    end
end












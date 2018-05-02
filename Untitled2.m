paws = 2:3;
trial = 2;

obsLimsY = [-0.0381 0.0381];
obsLimsZ = [0 .009];

close all; figure;
cmap = hsv(4);

for i = paws
    plot3(data(trial).locations(:,1,i), data(trial).locations(:,2,i), data(trial).locations(:,3,i), ...
        'color', cmap(i,:), 'linewidth', 3); hold on
%     plot(data(1).locations(:,1,i), data(1).locations(:,3,i), 'color', cmap(i,:)); hold on
end

patch(zeros(1,4), [obsLimsY(1) obsLimsY(1) obsLimsY(2) obsLimsY(2)], [obsLimsZ(1) obsLimsZ(2) obsLimsZ(2) obsLimsZ(1)], [0 0 0])
daspect([1 1 1]); pimpFig



close all;
figure('Position', [2164 355 560 420]);
p = plot(1:10, rand(1,10), 'linewidth', 5);
pause(.01)
set(p.Edge, 'ColorBinding', 'interpolated', 'ColorData', uint8([hsv(10)*255, ones(10,1)*120].'))
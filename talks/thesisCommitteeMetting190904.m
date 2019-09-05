%% global settings

%% formatting test

close all; figure('units', 'pixels', 'position', [100 100 500 300], 'color', 'black', 'menubar', 'none', 'inverthardcopy', false)
barFancy(rand(3,10), 'scatterAlpha', .8, 'colors', 'white')
print -clipboard -dbitmap

%%
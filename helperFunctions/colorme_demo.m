close all
numOffsets = 10;
colors = 2;

figHgt = 800/numOffsets;
offsets = linspace(0,1,numOffsets+1);
for i = 1:(length(offsets)-1)    
    colorme(colors, 'offset', offsets(i), 'showSamples', true, 'bgColor', 'white');
    set(gcf, 'position', [0 (i-1)*figHgt 400 figHgt])
end
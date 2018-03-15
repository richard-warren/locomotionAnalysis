

smps = 100;
x = linspace(0,10,1000);

y = repmat(x,smps,1)*.5 + randn(smps,length(x));
y2 = repmat(x,smps,1)*1 + randn(smps,length(x));

close all; figure;
shadedErrorBar(x, y, {@mean,@std}, 'lineProps', {'color', [0 0 1]}); hold on
shadedErrorBar(x, y2, {@mean,@std}, 'lineProps', {'color', [1 0 0]});
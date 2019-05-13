%% MAKE FAKE NONLINEAR DATA

x = 0:(2*pi/100):2*pi*5;
f = 1;
p = pi;
a = 1;
noise = rand(size(x)) * .0;

y1 = a*(sin(f*x) + f*x) + noise;
y2 = a*(sin(f*x+p) + f*x) + noise;
y = max(y1,y2);


% close all; figure; hold on;
% plot(x,y1);
% plot(x,y2);
% plot(x, y, 'LineWidth', 4, 'Color', [0 0 0 .4])
% pimpFig


% FIT MODEL

% paramsInit = rand(1,3);
paramsInit = [f p a] + rand(1,3);
fun = @(paramsInit) sseval(paramsInit, x, y);
paramsFitted = fminsearch(fun, paramsInit);
[~,yFit] = sseval(paramsFitted, x, y);

close all; figure; hold on
plot(x,y);
plot(x,yFit,'LineStyle', '--');
pimpFig


%% MAKE FAKE LINEAR DATA

x = 0:(2*pi/100):2*pi*5;
f = 1+rand(1);
p = pi+rand(1);
a = 1+rand(1);
noise = rand(size(x)) * .0;

y = a*(sin(f*x+p) + f*x) + noise;


% close all; figure; hold on;
% plot(x, y, 'LineWidth', 4, 'Color', [0 0 0 .4])
% pimpFig


% FIT MODEL
paramsInit = rand(1,3);
fun = @(f,p,a,x) a*(sin(f*x+p) + f*x);
funFit = fittype(fun, 'independent', {'x'});
paramsFitted = fit(x', y', funFit, ...
    'Lower', [-inf -pi -inf], 'Upper', [inf pi inf], ...
    'MaxIter', 10000);

close all;
figure; plot(paramsFitted,x,y); pimpFig








